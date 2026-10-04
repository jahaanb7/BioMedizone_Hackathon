"""s2: overlap variants with active plasma-cell/myeloma regulatory elements.

Adds one boolean column per registered peak track (`in_<track>`), the generic
cCRE baseline (`in_encode_ccre`), and a documented weighted
`regulatory_context_score`. Variants with no active-element overlap are
excluded with reason `no_active_overlap` (written to CSV).
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

# categories that count as "active plasma-cell regulatory elements"
ACTIVE_WEIGHT_KEYS = ("accessibility", "h3k27ac", "h3k4me1", "h3k4me3",
                      "super_enhancer", "tf_peak")


def step_regulatory_overlap(df: pd.DataFrame, *, bundle: ReferenceBundle,
                            config: AppConfig,
                            out_dir: Path) -> tuple[pd.DataFrame, FunnelStep, list[str]]:
    """Add track overlaps + regulatory_context_score; keep active-element hits."""
    rows_in = len(df)
    warnings: list[str] = []
    df = df.copy()
    df["start0"] = df["pos"].astype(np.int64) - 1

    tracks = bundle.peak_tracks()
    if not tracks:
        raise RuntimeError("no peak tracks registered (data/tracks/tracks.json); run make data")

    active_track_names: list[str] = []
    cat_tracks: dict[str, list[str]] = {k: [] for k in ACTIVE_WEIGHT_KEYS}
    active_lists = [[] for _ in range(len(df))]

    for t in tracks:
        idx = bundle.track_index(t["name"])
        col = f"in_{t['name']}"
        mask = np.zeros(len(df), dtype=bool)
        pos0 = df["start0"].to_numpy()
        for chrom, gi in df.groupby("chrom", sort=False).indices.items():
            mask[gi] = idx.contains_mask(chrom, pos0[gi])
        df[col] = mask
        cat = t.get("category", "other_histone")
        if cat in cat_tracks:
            cat_tracks[cat].append(col)
            active_track_names.append(t["name"])
        if cat == "super_enhancer":
            df["in_super_enhancer"] = df["in_super_enhancer"].fillna(False) | mask
        elif "in_super_enhancer" not in df.columns:
            df["in_super_enhancer"] = False
        if mask.any():
            for gi in np.nonzero(mask)[0]:
                active_lists[gi].append(t["name"])
    if "in_super_enhancer" not in df.columns:
        df["in_super_enhancer"] = False

    # generic cCRE baseline (always recorded, never used to drop variants)
    ccre = bundle.ccre()
    ccre_mask = np.zeros(len(df), dtype=bool)
    pos0 = df["start0"].to_numpy()
    for chrom, gi in df.groupby("chrom", sort=False).indices.items():
        ccre_mask[gi] = ccre.contains_mask(chrom, pos0[gi])
    df["in_encode_ccre"] = ccre_mask

    # weighted regulatory context score over categories that have data
    weights = config.regulatory_context_weights
    components: dict[str, np.ndarray] = {}
    for cat in ACTIVE_WEIGHT_KEYS:
        cols = cat_tracks[cat]
        if cat == "super_enhancer":
            components[cat] = df["in_super_enhancer"].fillna(False).to_numpy(bool)
        elif cols:
            m = np.zeros(len(df), dtype=bool)
            for c in cols:
                m |= df[c].to_numpy(bool)
            components[cat] = m
        elif cat in weights:
            warnings.append(f"category '{cat}' has no registered track; excluded from "
                            "regulatory_context_score denominator")
    available = {c: w for c, w in weights.items() if c in components}
    denom = sum(available.values()) or 1.0
    score = np.zeros(len(df))
    for c, w in available.items():
        score += w * components[c].astype(float)
    df["regulatory_context_score"] = score / denom
    df["active_tracks"] = [",".join(x) for x in active_lists]

    # keep variants overlapping >= 1 active element
    active_any = np.zeros(len(df), dtype=bool)
    for c in components:
        active_any |= components[c]
    drop = ~active_any
    n_drop = int(drop.sum())
    if n_drop:
        df.loc[drop].head(200_000).to_csv(out_dir / "s2_no_active_overlap.csv", index=False)
        df = df.loc[~drop].reset_index(drop=True)

    used = [{
        "name": t["name"], "display_name": t.get("display_name", t["name"]),
        "experiment": t.get("experiment", ""), "file_accession": t.get("file_accession", ""),
        "assay": t.get("assay", ""), "target": t.get("target", ""),
        "cell_type": t.get("cell_type", ""), "category": t.get("category", ""),
        "source_url": t.get("url", ""),
    } for t in tracks]
    (out_dir / "s2_tracks_used.json").write_text(json.dumps(used, indent=2))

    df = df.drop(columns=["start0"])
    df.to_parquet(out_dir / "s2_regulatory.parquet", index=False)
    step = FunnelStep(
        step="s2", label="Active regulatory element overlap", rows_in=rows_in,
        rows_out=len(df), removed=rows_in - len(df),
        removed_reasons={"no_active_overlap": n_drop} if n_drop else {},
        note=f"{len(tracks)} tracks; context score = "
             + ", ".join(f"{c}={w}" for c, w in available.items()))
    log.info("s2: %d -> %d across %d tracks", rows_in, len(df), len(tracks))
    return df, step, warnings
