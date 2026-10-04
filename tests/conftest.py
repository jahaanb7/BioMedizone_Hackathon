"""Shared fixtures: handcrafted PWM, tiny FASTA, fake ReferenceBundle for the motif step."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

CONSENSUS = "CACGTGAC"          # 8-mer used by the handcrafted PWM
PWM_W = len(CONSENSUS)


def make_pwm(consensus: str = CONSENSUS, match: float = 2.0,
             mismatch: float = -2.0) -> np.ndarray:
    """One-hot log-odds-style PWM, (4, W) ACGT-row layout as reference.py
    builds them: +match for the consensus base, else mismatch."""
    bases = {"A": 0, "C": 1, "G": 2, "T": 3}
    w = len(consensus)
    pwm = np.full((4, w), mismatch, dtype=np.float64)
    for i, ch in enumerate(consensus):
        pwm[bases[ch], i] = match
    return pwm


# Loci inside the tiny genome (0-based):
#   [100, 108) perfect consensus  -> SNV at 103 (G>T) destroys it  => broken
#   [250, 258) consensus with 1 mismatch (A at 253) -> A>G completes it => created
#   pos 300: df REF deliberately disagrees with the FASTA  => mismatch skip
#   [350, 358) revcomp(consensus) = opposite-strand site -> A>T at 353 => broken
TINY_GENOME_LEN = 400

RC_CONSENSUS = "GTCACGTG"       # revcomp of CONSENSUS


def build_tiny_fasta(path: Path) -> Path:
    seq = ["G"] * TINY_GENOME_LEN
    for i, ch in enumerate(CONSENSUS):
        seq[100 + i] = ch
    ref2 = "CACATGAC"                      # 1 mismatch vs consensus
    for i, ch in enumerate(ref2):
        seq[250 + i] = ch
    for i, ch in enumerate(RC_CONSENSUS):   # motif on the opposite strand
        seq[350 + i] = ch
    text = "".join(seq)
    fa = path / "tiny.fa"
    with open(fa, "w") as fh:
        fh.write(">chrTest\n")
        for i in range(0, len(text), 60):
            fh.write(text[i:i + 60] + "\n")
    pysam.faidx(str(fa))
    return fa


class _FakeBundle:
    """Duck-typed ReferenceBundle exposing only what step_motif_effect uses."""

    def __init__(self, motifs: dict, fasta_path: Path):
        self._motifs = motifs
        self._fa = pysam.FastaFile(str(fasta_path))

    def motifs(self) -> dict:
        return self._motifs

    def has_fasta(self, chrom: str) -> bool:
        return chrom in self._fa.references

    def fasta(self, chrom: str):
        return self._fa


@pytest.fixture()
def tiny_fasta(tmp_path: Path) -> Path:
    return build_tiny_fasta(tmp_path)


@pytest.fixture()
def motif_bundle(tiny_fasta: Path) -> _FakeBundle:
    pwm = make_pwm()
    # max attainable score = sum over positions of the best base
    motifs = {"TESTF": {"pwm": pwm,
                        "max_score": float(pwm.max(axis=0).sum()),
                        "matrix_id": "MA9999.1"}}
    return _FakeBundle(motifs, tiny_fasta)


@pytest.fixture()
def motif_frame() -> pd.DataFrame:
    """Four rows: broken, created, REF/FASTA mismatch, opposite-strand broken."""
    df = pd.DataFrame([
        # pos is 1-based VCF; variant sits inside the perfect site at 0-based 103
        {"variant_id": "chrTest:104:G:T", "chrom": "chrTest", "pos": 104,
         "ref": "G", "alt": "T"},
        # completes the 1-mismatch site at 0-based 253
        {"variant_id": "chrTest:254:A:G", "chrom": "chrTest", "pos": 254,
         "ref": "A", "alt": "G"},
        # REF disagrees with the tiny FASTA (0-based 300 holds 'G')
        {"variant_id": "chrTest:301:A:G", "chrom": "chrTest", "pos": 301,
         "ref": "A", "alt": "G"},
        # destroys the opposite-strand site at [350, 358) (A at 0-based 353)
        {"variant_id": "chrTest:354:A:T", "chrom": "chrTest", "pos": 354,
         "ref": "A", "alt": "T"},
        # 2 bp deletion removing motif bases at [101,103) -> mixed window
        # lengths (ref window 63, alt window 61) exercise length grouping
        {"variant_id": "chrTest:101:CAC:C", "chrom": "chrTest", "pos": 101,
         "ref": "CAC", "alt": "C"},
    ])
    return df.reset_index(drop=True)
