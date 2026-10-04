"""DepMap summaries: myeloma-vs-other CRISPR gene effect + expression flags.

Column names and lineage values are discovered at runtime from the downloaded
release files (they change between quarterly releases).
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

MM_PATTERN = re.compile(r"myeloma|plasma cell", re.I)


def find_depmap_files(cache_dir: Path | str) -> dict[str, Path]:
    """Locate Model / gene-effect / expression CSVs in the DepMap cache dir."""
    d = Path(cache_dir) / "depmap"
    if not d.exists():
        raise FileNotFoundError(f"DepMap cache {d} missing; run `make data`")
    out: dict[str, Path] = {}
    for p in sorted(d.glob("*.csv")):
        name = p.name.lower()
        if name == "model.csv":
            out["model"] = p
        elif "geneeffect" in name.replace("_", "").replace("-", ""):
            out["gene_effect"] = p
        elif "expression" in name:
            out["expression"] = p
        elif "somaticmutation" in name.replace("_", ""):
            out["somatic"] = p
    if "model" not in out:
        raise FileNotFoundError(f"Model.csv not found in {d}: {list(d.glob('*'))}")
    return out


def myeloma_model_ids(model_path: Path) -> tuple[list[str], str]:
    """(model_ids, description) for plasma-cell-myeloma lines in Model.csv."""
    df = pd.read_csv(model_path, low_memory=False)
    id_col = next((c for c in df.columns if c.lower() in ("modelid", "model_id")), None)
    if id_col is None:
        raise RuntimeError(f"Model.csv has no ModelID column: {list(df.columns)}")
    scan_cols = [c for c in df.columns
                 if re.search(r"oncotree|lineage|disease|subtype|primary", c, re.I)]
    if not scan_cols:
        raise RuntimeError(f"no lineage/subtype columns found in Model.csv: {list(df.columns)}")
    mask = pd.Series(False, index=df.index)
    hits: list[str] = []
    for c in scan_cols:
        m = df[c].astype(str).str.contains(MM_PATTERN, na=False)
        if m.any():
            mask |= m
            hits.append(f"{c}={sorted(df.loc[m, c].unique())[:3]}")
    ids = df.loc[mask, id_col].astype(str).tolist()
    if len(ids) < 3:
        raise RuntimeError(f"only {len(ids)} myeloma lines matched from columns {scan_cols}")
    desc = "; ".join(hits)
    log.info("DepMap myeloma lines: %d (%s)", len(ids), desc)
    return ids, desc


def _symbol(col: str) -> str:
    return col.split(" (")[0].removesuffix("_measured").strip()


def build_depmap_summary(cache_dir: Path | str, *, mm_effect_threshold: float = -0.5,
                         selectivity_margin: float = 0.3,
                         expressed_pct_threshold: int = 50,
                         force: bool = False) -> pd.DataFrame:
    """Per-gene: mean Chronos effect in myeloma vs other lines, selective and
    expressed flags. Cached to <cache_dir>/depmap/depmap_summary.parquet."""
    files = find_depmap_files(cache_dir)
    cache = files["model"].parent / "depmap_summary.parquet"
    if cache.exists() and not force:
        return pd.read_parquet(cache)

    mm_ids, mm_desc = myeloma_model_ids(files["model"])
    mm_set = set(mm_ids)

    ge = pd.read_csv(files["gene_effect"], index_col=0, low_memory=False)
    # DepMap releases vary: either models x genes or genes x models. Orient so
    # model IDs end up as COLUMNS (discovered at runtime from Model.csv IDs).
    if set(map(str, ge.index)) & mm_set:
        ge = ge.T
    mm_cols = [c for c in ge.columns if str(c) in mm_set]
    other_cols = [c for c in ge.columns if str(c) not in mm_set]
    if len(mm_cols) < 3:
        raise RuntimeError(f"gene-effect matrix has {len(mm_cols)} myeloma columns (need >=3); "
                           f"index sample: {list(ge.index[:3])} cols: {list(ge.columns[:3])}")
    mm_mean = ge[mm_cols].mean(axis=1, skipna=True)
    other_mean = ge[other_cols].mean(axis=1, skipna=True)
    summary = pd.DataFrame({
        "gene": [_symbol(c) for c in ge.index],
        "depmap_mm_mean_effect": mm_mean.to_numpy(),
        "depmap_other_mean_effect": other_mean.to_numpy(),
    })
    summary["depmap_mm_selective"] = (
        (summary["depmap_mm_mean_effect"] <= mm_effect_threshold)
        & ((summary["depmap_other_mean_effect"] - summary["depmap_mm_mean_effect"])
           >= selectivity_margin))
    log.info("DepMap gene effect: %d genes; %d myeloma-selective (thr=%.2f, margin=%.2f); "
             "myeloma lines=%d via {%s}", len(summary), int(summary["depmap_mm_selective"].sum()),
             mm_effect_threshold, selectivity_margin, len(mm_cols), mm_desc)

    if "expression" in files:
        ex = pd.read_csv(files["expression"], index_col=0, low_memory=False)
        if set(map(str, ex.index)) & mm_set:
            ex = ex.T
        ex_mm = [c for c in ex.columns if str(c) in mm_set]
        if ex_mm:
            expr_vals = ex[ex_mm]
            expr_vals.index = [_symbol(i) for i in expr_vals.index]
            frac = (expr_vals >= 1.0).mean(axis=1)
            summary["expressed_in_mm"] = summary["gene"].map(
                frac).ge(expressed_pct_threshold / 100.0).fillna(False)
            log.info("DepMap expression: %d/%d genes expressed (>=1 in >=%d%% mm lines)",
                     int(summary['expressed_in_mm'].sum()), len(summary), expressed_pct_threshold)
        else:
            summary["expressed_in_mm"] = False
            log.warning("no myeloma expression columns found")
    else:
        summary["expressed_in_mm"] = False
        log.warning("no expression file; expressed_in_mm=False for all genes")

    summary = summary.drop_duplicates(subset="gene")
    summary.to_parquet(cache, index=False)
    return summary


def somatic_mutation_lines(files: dict[str, Path]) -> pd.DataFrame | None:
    """Mode C helper: mutation table restricted to myeloma models (if downloaded)."""
    if "somatic" not in files:
        return None
    mm_ids, _ = myeloma_model_ids(files["model"])
    mm_set = set(mm_ids)
    mut = pd.read_csv(files["somatic"], low_memory=False)
    id_col = next((c for c in mut.columns if re.search(r"modelid|model_id", c, re.I)), None)
    if id_col is None:
        raise RuntimeError(f"no ModelID column in somatic file: {list(mut.columns)}")
    sub = mut[mut[id_col].astype(str).isin(mm_set)]
    log.info("Mode C: %d/%d somatic mutations in %d myeloma lines", len(sub), len(mut), len(mm_set))
    return sub
