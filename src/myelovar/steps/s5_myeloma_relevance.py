"""s5: myeloma relevance scoring from DepMap, known-driver lists, expression
and GWAS-locus membership. Missing sources are recorded, never invented."""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from myelovar.config import AppConfig
from myelovar.reference import ReferenceBundle
from myelovar.schemas import FunnelStep

log = logging.getLogger(__name__)


def step_myeloma_relevance(df: pd.DataFrame, *, bundle: ReferenceBundle,
                           config: AppConfig,
                           out_dir: Path) -> tuple[pd.DataFrame, FunnelStep, list[str]]:
    """Attach DepMap / known-gene / expression / GWAS-locus columns + score."""
    rows_in = len(df)
    warnings: list[str] = []
    df = df.copy()

    # primary gene = first target gene, else nearest gene
    primary = df["target_genes"].fillna("").str.split(",").str[0].replace("", np.nan)
    primary = primary.fillna(df["nearest_gene"].replace("", np.nan))

    depmap = bundle.depmap()
    dm = depmap.set_index("gene")
    df["depmap_mm_mean_effect"] = primary.map(dm["depmap_mm_mean_effect"])
    df["depmap_mm_selective"] = primary.map(dm["depmap_mm_selective"]).fillna(False)
    df["expressed_in_mm"] = primary.map(dm["expressed_in_mm"]).fillna(False)

    mm_genes = bundle.mm_genes()
    known_set = set(mm_genes["gene"].astype(str)) if len(mm_genes) else set()
    src_map = (dict(zip(mm_genes["gene"].astype(str), mm_genes["source"].astype(str)))
               if len(mm_genes) and "source" in mm_genes.columns else {})
    df["is_known_mm_gene"] = primary.isin(known_set)
    df["known_gene_source"] = primary.map(src_map).fillna("")
    if not known_set:
        warnings.append("known-myeloma-gene list unavailable; is_known_mm_gene=False "
                        "(see data manifest for the failure)")

    # GWAS locus membership (built from Mode B / cached loci)
    loci = bundle.gwas_loci()
    in_locus = np.zeros(len(df), dtype=bool)
    if loci:
        loc_by_chrom: dict[str, list[tuple[int, int]]] = {}
        for L in loci:
            loc_by_chrom.setdefault(L["chrom"], []).append((L["start"], L["end"]))
        for chrom, gi in df.groupby("chrom", sort=False).indices.items():
            spans = loc_by_chrom.get(chrom)
            if not spans:
                continue
            starts = np.array([s for s, _ in spans])
            ends = np.array([e for _, e in spans])
            order = np.argsort(starts)
            starts, ends = starts[order], ends[order]
            pos0 = df.iloc[gi]["pos"].to_numpy(np.int64) - 1
            idx = np.searchsorted(starts, pos0, side="right") - 1
            valid = idx >= 0
            m = np.zeros(len(pos0), dtype=bool)
            m[valid] = pos0[valid] < ends[idx[valid]]
            in_locus[gi] = m
    else:
        warnings.append("no GWAS loci cached (gwas/loci.json); in_mm_gwas_locus=False for all")
    df["in_mm_gwas_locus"] = in_locus
    if "gwas_flag" not in df.columns:
        df["gwas_flag"] = False
    df["gwas_flag"] = df["gwas_flag"].fillna(False).astype(bool)

    # ---- weighted relevance score over *available* components ------------
    w = config.myeloma_relevance_weights
    comp: dict[str, np.ndarray] = {}
    depmap_present = df["depmap_mm_mean_effect"].notna()
    if depmap_present.any():
        dep_score = np.clip(-df["depmap_mm_mean_effect"].fillna(0.0), 0, 1.5) / 1.5
        dep_score = np.where(df["depmap_mm_selective"], np.maximum(dep_score, 0.7), dep_score)
        comp["depmap"] = np.asarray(dep_score, dtype=float) * \
            depmap_present.to_numpy(dtype=float)
    else:
        warnings.append("DepMap gene-effect unavailable for all variants; depmap component dropped")
    if known_set:
        comp["known_gene"] = df["is_known_mm_gene"].astype(float).to_numpy()
    if df["expressed_in_mm"].any():
        comp["expressed"] = df["expressed_in_mm"].astype(float).to_numpy()
    if loci:
        comp["gwas"] = (df["in_mm_gwas_locus"] | df["gwas_flag"]).astype(float).to_numpy()

    available_w = {k: w[k] for k in comp}
    denom = sum(available_w.values()) or 1.0
    score = np.zeros(len(df))
    for k, wk in available_w.items():
        score += wk * comp[k]
    df["myeloma_relevance_score"] = score / denom

    df.to_parquet(out_dir / "s5_relevance.parquet", index=False)
    step = FunnelStep(
        step="s5", label="Myeloma relevance", rows_in=rows_in, rows_out=len(df), removed=0,
        note="components: " + ", ".join(f"{k}={v}" for k, v in available_w.items())
             + (f"; {len(loci)} GWAS loci" if loci else ""))
    log.info("s5: %d rows; components=%s", rows_in, list(available_w))
    return df, step, warnings
