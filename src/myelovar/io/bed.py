"""Interval arithmetic: BED loading, merging, boolean/range queries.

All coordinates flowing through the pipeline are 0-based half-open [start, end)
internally; VCF 1-based positions are converted explicitly at the boundary
(see vcf_to_0based) and unit-tested.
"""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

OPENERS = {".gz": gzip.open, ".bgz": gzip.open}


def open_text(path: Path | str):
    p = Path(path)
    opener = OPENERS.get(p.suffix, open)
    return opener(p, "rt")


def vcf_to_0based(pos1: int, ref: str) -> tuple[int, int]:
    """Convert a VCF 1-based POS + REF allele to 0-based half-open interval."""
    return pos1 - 1, pos1 - 1 + len(ref)


def load_bed(path: Path | str, min_cols: int = 3) -> pd.DataFrame:
    """Load a BED-like file (comment lines skipped) into chrom/start/end frame."""
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if not line or line[0] == "#" or line.startswith("track") or line.startswith("browser"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < min_cols:
                continue
            rows.append((f[0], int(f[1]), int(f[2]), f[3:] or ()))
    if not rows:
        raise ValueError(f"no intervals loaded from {path}")
    df = pd.DataFrame([(c, s, e) for c, s, e, _ in rows], columns=["chrom", "start", "end"])
    extra = [r[3] for r in rows]
    return df, extra


def merge_intervals(df: pd.DataFrame) -> pd.DataFrame:
    """Merge overlapping intervals per chromosome (exact for boolean overlap)."""
    out = []
    for chrom, sub in df.groupby("chrom", sort=False):
        sub = sub.sort_values("start")
        starts = sub["start"].to_numpy()
        ends = sub["end"].to_numpy()
        cs, ce = starts[0], ends[0]
        for s, e in zip(starts[1:], ends[1:]):
            if s <= ce:
                ce = max(ce, e)
            else:
                out.append((chrom, cs, ce))
                cs, ce = s, e
        out.append((chrom, cs, ce))
    return pd.DataFrame(out, columns=["chrom", "start", "end"])


class IntervalIndex:
    """Per-chromosome sorted interval index with O(log n) membership queries."""

    def __init__(self, df: pd.DataFrame, merge: bool = True):
        if merge:
            df = merge_intervals(df)
        self.chroms: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for chrom, sub in df.groupby("chrom", sort=False):
            sub = sub.sort_values("start")
            self.chroms[chrom] = (sub["start"].to_numpy(np.int64),
                                  sub["end"].to_numpy(np.int64))

    def __len__(self) -> int:
        return int(sum(len(v[0]) for v in self.chroms.values()))

    def contains_mask(self, chrom: str, pos0: np.ndarray) -> np.ndarray:
        """Boolean mask: does 0-based position fall in any interval?"""
        out = np.zeros(len(pos0), dtype=bool)
        arr = self.chroms.get(chrom)
        if arr is None or len(pos0) == 0:
            return out
        starts, ends = arr
        idx = np.searchsorted(starts, pos0, side="right") - 1
        valid = idx >= 0
        out[valid] = pos0[valid] < ends[idx[valid]]
        return out

    def overlapping_intervals(self, chrom: str, start0: int, end0: int) -> list[tuple[int, int]]:
        """All intervals overlapping [start0, end0) (bounded back-scan)."""
        arr = self.chroms.get(chrom)
        if arr is None:
            return []
        starts, ends = arr
        i = int(np.searchsorted(starts, end0, side="left")) - 1
        hits: list[tuple[int, int]] = []
        while i >= 0 and starts[i] < end0 and len(hits) < 10_000:
            if ends[i] > start0:
                hits.append((int(starts[i]), int(ends[i])))
            i -= 1
        return hits

    def max_reach(self) -> int:
        """Max interval length (bounds exact back-scan depth)."""
        return int(max((e - s).max() for s, e in self.chroms.values()
                       if len(s) > 0), default=0)


def point_query_mask(index: IntervalIndex, chrom: str, pos0: np.ndarray) -> np.ndarray:
    return index.contains_mask(chrom, pos0)
