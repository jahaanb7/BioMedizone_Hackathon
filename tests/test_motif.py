"""Motif-effect step: handcrafted sequences with a known perfect/disrupted site.

Covers the low-level scorer (vs. an independent brute-force reference), the
reverse-strand variant-index math, and the full step_motif_effect on a tiny
FASTA where the expected outcome of each variant is known analytically.
"""
from __future__ import annotations

import numpy as np
import pytest

from myelovar.config import AppConfig
from myelovar.steps.s3_motif_effect import _best_scores, _encode, step_motif_effect

from conftest import CONSENSUS, make_pwm

BASES = {"A": 0, "C": 1, "G": 2, "T": 3}


# ------------------------------------------------------------------ _encode
def test_encode_maps_bases_and_masks_ambiguous():
    arr = _encode(["ACGTN"])
    assert arr.tolist() == [[0, 1, 2, 3, -1]]
    assert _encode([]).shape == (0, 0)


# ----------------------------------------------------------- _best_scores --
def _brute_best(seq: str, pwm: np.ndarray, varidx: int) -> float:
    """Independent reference: best score over sites containing varidx only."""
    w = pwm.shape[1]
    best = -np.inf
    for s in range(0, len(seq) - w + 1):
        if not (s <= varidx < s + w):
            continue
        sc = 0.0
        for i in range(w):
            b = BASES.get(seq[s + i])
            sc += 0.0 if b is None else pwm[b, i]
        best = max(best, sc)
    return best


def test_best_scores_matches_brute_force_randomized():
    rng = np.random.default_rng(11)
    pwm = make_pwm()
    w = pwm.shape[1]
    for _ in range(25):
        L = int(rng.integers(w, 60))
        seq = "".join(rng.choice(list("ACGTN"), size=L))
        varidx = int(rng.integers(0, L))
        got = _best_scores(_encode([seq]), pwm, np.array([varidx]))[0]
        want = _brute_best(seq, pwm, varidx)
        assert np.isclose(got, want), (seq, varidx, got, want)


def test_best_scores_ignores_perfect_site_outside_variant_window():
    # perfect site at [0,8); variant far to its right -> site must be excluded
    seq = CONSENSUS + "AAAAAA"          # len 14, variant at index 13
    pwm = make_pwm()
    got = _best_scores(_encode([seq]), pwm, np.array([13]))[0]
    perfect = float(pwm.max(axis=0).sum())
    assert got < perfect
    # with the variant inside the site, the perfect score is found
    got_in = _best_scores(_encode([seq]), pwm, np.array([3]))[0]
    assert np.isclose(got_in, perfect)


def test_reverse_strand_motif_detected_with_rc_pwm():
    """A site on the opposite strand (revcomp of the consensus) must be found
    by scoring the FORWARD sequence with the reverse-complement PWM - the
    exact scan step_motif_effect performs."""
    pwm = make_pwm()
    W = 30
    rc_consensus = "".join({"A": "T", "C": "G", "G": "C", "T": "A"}[b]
                           for b in reversed(CONSENSUS))
    # window length 2*W+1 (SNV REF length 1); opposite-strand motif covers the
    # variant at index W (motif start = W-3)
    seq = ["G"] * (2 * W + 1)
    start = W - 3
    for i, ch in enumerate(rc_consensus):
        seq[start + i] = ch
    seq = "".join(seq)
    assert seq[W] == rc_consensus[3]           # variant base inside the site
    varidx = np.array([W])

    got = _best_scores(_encode([seq]), pwm[::-1, ::-1], varidx)[0]
    perfect = float(pwm.max(axis=0).sum())
    assert np.isclose(got, perfect), "motif on the reverse strand was missed"
    # ...and the forward PWM does not score this site perfectly (specificity)
    fwd = _best_scores(_encode([seq]), pwm, varidx)[0]
    assert fwd < perfect


# ------------------------------------------------------ full step (FASTA) --
def test_step_motif_effect_known_outcomes(tmp_path, motif_bundle, motif_frame):
    cfg = AppConfig()
    out = tmp_path / "work"
    out.mkdir()
    df, step, warnings = step_motif_effect(motif_frame, bundle=motif_bundle,
                                           config=cfg, out_dir=out)
    # never drops rows
    assert len(df) == len(motif_frame)
    assert step.rows_in == step.rows_out == 5 and step.removed == 0
    assert (out / "s3_motif.parquet").exists()

    broken = df.iloc[0]
    assert broken["motif_effect"] == "broken"
    assert broken["top_motif_tf"] == "TESTF"
    assert np.isclose(broken["motif_ref_rel"], 1.0)       # perfect ref site
    assert np.isclose(broken["motif_alt_rel"], 0.75)       # 1 mismatch: (16-4)/16
    assert np.isclose(broken["motif_delta"], -0.25)
    assert np.isclose(broken["motif_disruption_score"], 0.25)

    created = df.iloc[1]
    assert created["motif_effect"] == "created"
    assert np.isclose(created["motif_ref_rel"], 0.75)
    assert np.isclose(created["motif_alt_rel"], 1.0)
    assert np.isclose(created["motif_delta"], 0.25)

    mismatch = df.iloc[2]
    assert mismatch["motif_effect"] == "none"
    assert mismatch["top_motif_tf"] is None or pd_isna(mismatch["top_motif_tf"])
    assert bool(mismatch["motif_scanned"]) is False
    assert len(warnings) == 1 and "mismatch" in warnings[0]

    rc_broken = df.iloc[3]                    # opposite-strand site destroyed
    assert rc_broken["motif_effect"] == "broken"
    assert rc_broken["top_motif_tf"] == "TESTF"
    assert np.isclose(rc_broken["motif_ref_rel"], 1.0)
    assert np.isclose(rc_broken["motif_alt_rel"], 0.75)
    assert np.isclose(rc_broken["motif_delta"], -0.25)

    deletion = df.iloc[4]                     # mixed window lengths (63/61)
    assert deletion["motif_effect"] == "broken"
    assert np.isclose(deletion["motif_ref_rel"], 1.0)
    # brute force confirms the best variant-overlapping alt site scores 8/16
    assert np.isclose(deletion["motif_alt_rel"], 0.5)
    assert np.isclose(deletion["motif_delta"], -0.5)


def pd_isna(v) -> bool:
    import pandas as pd
    return bool(pd.isna(v))


def test_step_motif_effect_raises_without_fasta(tmp_path, motif_bundle):
    import pandas as pd
    cfg = AppConfig()
    out = tmp_path / "work2"
    out.mkdir()
    df = pd.DataFrame([{"variant_id": "chrZ:10:A:T", "chrom": "chrZ",
                        "pos": 10, "ref": "A", "alt": "T"}])
    with pytest.raises(FileNotFoundError, match="FASTA missing"):
        step_motif_effect(df, bundle=motif_bundle, config=cfg, out_dir=out)
