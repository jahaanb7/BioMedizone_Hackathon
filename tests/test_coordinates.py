"""Coordinate conversion and interval-index arithmetic (0-based half-open rules)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from myelovar.io.bed import IntervalIndex, merge_intervals, vcf_to_0based


def test_vcf_to_0based_snv():
    # VCF POS is 1-based; an SNV at POS=100 covers exactly base 99 (0-based)
    assert vcf_to_0based(100, "A") == (99, 100)


def test_vcf_to_0based_deletion_and_insertion():
    # deletion of 3 bases at POS=50 -> [49, 52)
    assert vcf_to_0based(50, "ACG") == (49, 52)
    # insertion: REF is a single anchor base
    start, end = vcf_to_0based(1, "A")
    assert (start, end) == (0, 1)


def test_vcf_to_0based_matches_single_base_interval():
    for pos in (1, 2, 1_000, 248_956_422):
        s, e = vcf_to_0based(pos, "T")
        assert e - s == 1 and s == pos - 1


def test_merge_intervals_overlapping_and_adjacent():
    df = pd.DataFrame([
        ("chr1", 10, 20), ("chr1", 15, 25),   # overlaps -> merge
        ("chr1", 30, 40), ("chr1", 40, 45),   # touching -> merge (half-open)
        ("chr1", 100, 110),                   # separate
        ("chr2", 1, 5),
    ], columns=["chrom", "start", "end"])
    m = merge_intervals(df)
    chr1 = m[m["chrom"] == "chr1"].values.tolist()
    assert chr1 == [["chr1", 10, 25], ["chr1", 30, 45], ["chr1", 100, 110]]
    assert len(m[m["chrom"] == "chr2"]) == 1


def _brute_contains(intervals, pos0):
    return any(s <= pos0 < e for s, e in intervals)


def test_interval_index_contains_mask_matches_brute_force():
    rng = np.random.default_rng(7)
    starts = np.sort(rng.integers(0, 1000, size=50))
    ends = starts + rng.integers(1, 60, size=50)
    df = pd.DataFrame({"chrom": ["chr7"] * 50, "start": starts, "end": ends})
    idx = IntervalIndex(df)  # merges overlapping intervals first
    merged = merge_intervals(df)
    pairs = list(zip(merged["start"], merged["end"]))

    probe = rng.integers(-5, 1010, size=200)
    got = idx.contains_mask("chr7", probe)
    want = np.array([_brute_contains(pairs, int(p)) for p in probe])
    assert np.array_equal(got, want)
    # unknown chromosome -> all False, no crash
    assert not idx.contains_mask("chrZ", probe).any()


def test_interval_index_overlapping_intervals_half_open_semantics():
    df = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [100, 200],
                       "end": [150, 250]})
    idx = IntervalIndex(df, merge=False)
    # query touching only the left edge of [100,150): [150,160) must NOT hit
    assert idx.overlapping_intervals("chr1", 150, 160) == []
    # [90,100) touches start -> no overlap (half-open)
    assert idx.overlapping_intervals("chr1", 90, 100) == []
    # [90,101) overlaps by one base -> hit
    assert idx.overlapping_intervals("chr1", 90, 101) == [(100, 150)]
