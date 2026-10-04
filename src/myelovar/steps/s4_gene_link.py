"""s4: link each variant to likely target genes.

Evidence types (confidence per spec):
  * ABC predictions (Nasser 2021, chosen biosample; proxy labelled in metadata)
  * distance to TSS within +/-100 kb with distance decay
  * nearest gene as last-resort fallback
  high = ABC + distance agree; medium = ABC only; low = distance only/nearest.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from myelovar.config import AppConfig
from myelovar.reference import ReferenceBundle
from myelovar.schemas import FunnelStep

log = logging.getLogger(__name__)


class _PayloadIndex:
    """Raw (unmerged) interval index returning payload row ids on overlap."""

    def __init__(self, df: pd.DataFrame, payload_col: str):
        self.chroms: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        self.max_width = 1
        for chrom, sub in df.groupby("chrom", sort=False):
            sub = sub.sort_values("start")
            self.chroms[chrom] = (sub["start"].to_numpy(np.int64),
                                  sub["end"].to_numpy(np.int64),
                                  sub[payload_col].to_numpy())
            self.max_width = max(self.max_width,
                                 int((sub["end"] - sub["start"]).max()))

    def overlap(self, chrom: str, start0: int, end0: int) -> np.ndarray:
        arr = self.chroms.get(chrom)
        if arr is None:
            return np.empty(0, dtype=np.int64)
        starts, ends, payload = arr
        i = int(np.searchsorted(starts, end0, side="left")) - 1
        hits: list[int] = []
        while i >= 0 and starts[i] > start0 - self.max_width and len(hits) < 500:
            if ends[i] > start0:
                hits.append(int(payload[i]))
                if starts[i] <= start0 and ends[i] > start0:
                    # exhaustive-enough: continue a bit for nested intervals
                    pass
            i -= 1
            if i >= 0 and starts[i] < start0 - self.max_width:
                break
        return np.asarray(hits, dtype=np.int64)


def step_gene_link(df: pd.DataFrame, *, bundle: ReferenceBundle,
                   config: AppConfig,
                   out_dir: Path) -> tuple[pd.DataFrame, FunnelStep, list[str]]:
    """Attach target_genes / link_* columns to every surviving variant."""
    rows_in = len(df)
    warnings: list[str] = []
    abc, abc_meta = bundle.abc()
    if abc_meta.get("proxy"):
        warnings.append(
            f"ABC predictions use biosample {abc_meta.get('cell_types')} as a labelled "
            "proxy (no plasma-cell biosample in Nasser et al. 2021)")

    abc_idx = _PayloadIndex(abc.assign(_row=np.arange(len(abc))), "_row")
    dist_kb = config.link.distance_kb
    scale_kb = config.link.decay_scale_kb

    target_genes: list[str] = []
    details_col: list[str] = []
    methods: list[str] = []
    ev_counts: list[int] = []
    confidences: list[str] = []
    link_scores: list[float] = []

    n_abc = n_dist = n_nearest = 0
    for row in df.itertuples(index=False):
        pos0 = int(row.pos) - 1
        chrom = row.chrom
        # --- ABC evidence
        abc_hits = abc_idx.overlap(chrom, pos0, pos0 + 1)
        abc_genes: dict[str, float] = {}
        if len(abc_hits):
            sub = abc.iloc[abc_hits]
            for g, s in zip(sub["gene"], sub["abc_score"]):
                g = str(g)
                abc_genes[g] = max(abc_genes.get(g, 0.0), float(s))
        # --- distance evidence (TSS within window)
        near = bundle.genes_tss_within(chrom, pos0, dist_kb * 1000)
        dist_genes: dict[str, tuple[int, float]] = {}
        for g, tss in zip(near["gene"], near["tss"]):
            d = abs(int(tss) - pos0)
            decay = 0.5 / (1.0 + (d / (scale_kb * 1000)) ** 2)
            key = str(g)
            if key not in dist_genes or d < dist_genes[key][0]:
                dist_genes[key] = (d, decay)

        # --- combine evidence per gene
        combined: list[dict] = []
        for g in set(abc_genes) | set(dist_genes):
            ev = []
            score = 0.0
            if g in abc_genes:
                ev.append("ABC")
                score = max(score, min(1.0, abc_genes[g]))
            detail = {"gene": g, "evidence": ev, "score": round(score, 4)}
            if g in dist_genes:
                d, decay = dist_genes[g]
                ev.append("distance")
                detail["tss_distance_bp"] = d
                detail["distance_score"] = round(decay, 4)
                score = max(score, decay)
                detail["evidence"] = ev
                detail["score"] = round(score, 4)
            if ev:
                combined.append(detail)
        combined.sort(key=lambda x: -x["score"])

        if combined:
            primary = combined[0]
            genes = [c["gene"] for c in combined[:5]]
            ev = primary["evidence"]
            if "ABC" in ev and "distance" in ev:
                conf = "high"
            elif "ABC" in ev:
                conf = "medium"
            else:
                conf = "low"
            method = "+".join(ev) if ev else "none"
            n_abc += int("ABC" in ev)
            n_dist += 1
            score = float(primary["score"])
        else:
            ng = row.nearest_gene if getattr(row, "nearest_gene", None) else ""
            genes = [ng] if ng else []
            primary = {"gene": ng, "evidence": ["nearest_gene"], "score": 0.0}
            conf, method, ev = "low", "nearest_gene", ["nearest_gene"]
            n_nearest += 1
            score = 0.0

        target_genes.append(",".join(g for g in genes if g))
        details_col.append(json.dumps(combined[:5]))
        methods.append(method)
        ev_counts.append(len(primary.get("evidence", [])))
        confidences.append(conf)
        link_scores.append(score)

    df = df.copy()
    df["target_genes"] = target_genes
    df["link_method"] = methods
    df["link_evidence_count"] = ev_counts
    df["link_confidence"] = confidences
    df["link_score"] = link_scores
    df["link_details"] = details_col
    df.to_parquet(out_dir / "s4_linked.parquet", index=False)
    step = FunnelStep(
        step="s4", label="Enhancer-to-gene linking", rows_in=rows_in, rows_out=len(df),
        removed=0,
        note=f"ABC biosample={abc_meta.get('cell_types')} proxy={abc_meta.get('proxy')}; "
             f"primary links: ABC={n_abc}, distance={n_dist}, nearest-fallback={n_nearest}")
    log.info("s4: %d variants linked (ABC=%d distance=%d nearest=%d)",
             rows_in, n_abc, n_dist, n_nearest)
    return df, step, warnings
