"""s1: quality -> blacklist -> common-variant -> coding filters + annotation.

Order is cheapest-first; every removal writes a reason-coded CSV; counts are
reported in the returned FunnelStep (fail loudly: raises on missing reference).

Adaptations for the GIAB benchmark VCF (documented per spec): its FORMAT is
GT:PS:DP:ADALL:AD:GQ, so DP>=10 and GQ>=20 apply directly; all records in the
first 200k were FILTER=PASS (high-confidence benchmark calls).
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from myelovar.config import AppConfig
from myelovar.io.bed import IntervalIndex, vcf_to_0based
from myelovar.reference import ReferenceBundle
from myelovar.schemas import FunnelStep

log = logging.getLogger(__name__)


def _per_chrom_mask(df: pd.DataFrame, fn) -> np.ndarray:
    """Apply fn(chrom, pos0_array) -> mask over rows of df (grouped by chrom)."""
    out = np.zeros(len(df), dtype=bool)
    pos0 = df["pos"].to_numpy(np.int64) - 1
    for chrom, idx in df.groupby("chrom", sort=False).indices.items():
        out[idx] = fn(chrom, pos0[idx])
    return out


def _nearest_gene_vectorised(df: pd.DataFrame, genes: pd.DataFrame,
                             max_back: int = 3) -> tuple[list[str], list[int]]:
    """Nearest gene by span using a bounded back-scan + next-gene candidate.

    Vectorised: for each variant take candidates {i, i-1..i-max_back, i+1}
    where i = last gene starting at/before the position; distance 0 if inside.
    """
    nearest = np.full(len(df), None, dtype=object)
    dists = np.full(len(df), -1, dtype=np.int64)
    pos0 = df["pos"].to_numpy(np.int64) - 1
    for chrom, idx in df.groupby("chrom", sort=False).indices.items():
        g = genes[genes["chrom"] == chrom]
        if g.empty:
            continue
        g = g.sort_values("start")
        starts = g["start"].to_numpy(np.int64)
        ends = g["end"].to_numpy(np.int64)
        names = g["gene"].to_numpy()
        p = pos0[idx]
        i = np.searchsorted(starts, p, side="right") - 1
        best_d = np.full(len(p), np.iinfo(np.int64).max, dtype=np.int64)
        best_g = np.full(len(p), None, dtype=object)

        def consider(c: np.ndarray) -> None:
            nonlocal best_d, best_g
            cc = np.clip(c, 0, len(starts) - 1)
            s, e = starts[cc], ends[cc]
            d = np.where(p < s, s - p, np.where(p >= e, p - e, 0))
            better = d < best_d
            best_d = np.where(better, d, best_d)
            best_g = np.where(better, names[cc], best_g)

        consider(i)
        for k in range(1, max_back + 1):
            consider(i - k)
        consider(i + 1)
        nearest[idx] = best_g
        dists[idx] = best_d
    return [g if g is not None else "" for g in nearest], dists.tolist()


def step_annotate_filter(df: pd.DataFrame, *, bundle: ReferenceBundle,
                         config: AppConfig, out_dir: Path,
                         mode: str = "A") -> tuple[pd.DataFrame, FunnelStep]:
    """Annotate and filter variants; returns (kept_df, FunnelStep)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    reasons: dict[str, int] = {}
    warnings: list[str] = []
    rows_in = len(df)
    df = df.copy()
    df["start0"], df["end0"] = vcf_to_0based_arr(df["pos"], df["ref"])

    # --- 1. quality ------------------------------------------------------
    if mode == "A":
        if "filter" in df.columns:
            ok = df["filter"].isin(["PASS", "."])
            _cull(df, ok, "not_PASS", out_dir / "s1_removed_quality.csv", reasons)
        for col, thr in (("dp", config.filter.min_dp), ("gq", config.filter.min_gq)):
            if col in df.columns and df[col].notna().any():
                ok = df[col].isna() | (df[col] >= thr)
                _cull(df, ok, f"{col}_lt_{thr}", out_dir / "s1_removed_quality.csv", reasons)
            else:
                warnings.append(f"input has no usable {col.upper()} values; "
                                f"{col.upper()}>={thr} filter skipped")
    else:
        warnings.append("Mode B: quality thresholds (DP/GQ/FILTER) not applicable "
                        "to API-derived locus variants")

    # --- 2. variant size (SNVs + small indels only) ----------------------
    ok = (df["end0"] - df["start0"]) <= 50
    _cull(df, ok, "variant_gt_50bp", out_dir / "s1_removed_large.csv", reasons)

    # --- 3. blacklist ----------------------------------------------------
    bl = bundle.blacklist()
    ok = ~_per_chrom_mask(df, lambda c, p: bl.contains_mask(c, p))
    _cull(df, ok, "blacklisted_region", out_dir / "s1_removed_blacklist.csv", reasons)

    # --- 4. common variants (Mode A only) --------------------------------
    af_cols = bundle.af(df)
    df = df.merge(af_cols, on="variant_id", how="left", validate="one_to_one")
    if mode == "A" and config.filter.remove_common:
        high = df["af"].notna() & (df["af"] > config.filter.max_af)
        track_only = df["af"].isna() & df["in_common_track"].fillna(False)
        ok = ~(high | track_only)
        df["_af_reason"] = np.where(high, f"af_gt_{config.filter.max_af}",
                                    np.where(track_only, "in_dbSNP_common_track", ""))
        _cull(df, ok, "common_variant", out_dir / "s1_removed_common.csv", reasons)
        df = df.drop(columns=["_af_reason"])
    elif mode == "A":
        warnings.append("common-variant filter disabled in config")
    else:
        warnings.append("Mode B: common-variant filter disabled (GWAS signals are common)")

    # --- 5. coding split + splice flag -----------------------------------
    cds = bundle.cds_index()
    exons = bundle.exon_index()
    cds_expanded = IntervalIndex(
        _expand(bundle.gencode()["cds"], 2), merge=True)
    in_cds = _per_chrom_mask(df, lambda c, p: cds.contains_mask(c, p))
    near_splice = _per_chrom_mask(df, lambda c, p: cds_expanded.contains_mask(c, p))
    in_exon = _per_chrom_mask(df, lambda c, p: exons.contains_mask(c, p))
    df["is_splice"] = near_splice & ~in_cds
    coding = in_cds & ~df["is_splice"].to_numpy()
    df["consequence_class"] = np.select(
        [coding, df["is_splice"], in_exon],
        ["coding", "splice", "exonic_noncoding"], default="noncoding")

    if coding.any():
        df[coding].drop(columns=["start0", "end0"], errors="ignore").to_csv(
            out_dir / "coding_variants.csv", index=False)
    _cull(df, ~coding, "coding_variant", out_dir / "s1_removed_coding.csv", reasons)

    # --- 6. nearest gene annotation --------------------------------------
    nearest, ndist = _nearest_gene_vectorised(df, bundle.genes())
    df["nearest_gene"] = nearest
    df["nearest_gene_dist"] = ndist

    df = df.drop(columns=["start0", "end0"], errors="ignore").reset_index(drop=True)
    df.to_parquet(out_dir / "s1_annotated.parquet", index=False)
    step = FunnelStep(step="s1", label="Annotate + filter (quality, blacklist, common, coding)",
                      rows_in=rows_in, rows_out=len(df),
                      removed=rows_in - len(df), removed_reasons=reasons,
                      note="; ".join(warnings) if warnings else "")
    log.info("s1: %d -> %d (removed %s)", rows_in, len(df), reasons)
    return df, step, warnings


def _cull(df: pd.DataFrame, keep_mask: np.ndarray, reason: str,
          path: Path, reasons: dict[str, int]) -> None:
    """Drop rows where keep_mask is False, appending them to a removal CSV."""
    drop = ~keep_mask
    n = int(drop.sum())
    if n:
        removed = df.loc[drop]
        header = not path.exists()
        removed.head(200_000).to_csv(path, mode="a", header=header, index=False)
        reasons[reason] = reasons.get(reason, 0) + n
        df.drop(df.index[drop], inplace=True)


def vcf_to_0based_arr(pos: pd.Series, ref: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Vectorised 1-based POS + REF -> 0-based half-open interval."""
    start = pos.astype(np.int64) - 1
    return start, start + ref.str.len().astype(np.int64)


def _expand(df: pd.DataFrame, by: int) -> pd.DataFrame:
    out = df.copy()
    out["start"] = (out["start"] - by).clip(lower=0)
    out["end"] = out["end"] + by
    return out
