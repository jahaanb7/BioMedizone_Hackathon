"""s3: lightweight TF-motif disruption scan (log-odds over JASPAR PWMs).

For every surviving variant we extract a +/-30 bp reference window (bgzip
FASTA via pysam), build the ALT window, and scan both strands of both alleles
for every prioritised TF matrix, considering only sites that *contain the
variant position*. Relative scores are vs. the PWM's maximum attainable score.

Flags (thresholds in config.motif):
  broken:  ref_rel >= min_strong and (ref_rel - alt_rel) >= min_drop
  created: alt_rel >= min_strong and (alt_rel - ref_rel) >= min_drop
Reference/FASTA mismatches are counted and skipped (build sanity check).
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from myelovar.config import AppConfig
from myelovar.reference import ReferenceBundle
from myelovar.schemas import FunnelStep

log = logging.getLogger(__name__)

BASES = {"A": 0, "C": 1, "G": 2, "T": 3}
FLANK = 30


def _encode(seqs: list[str]) -> np.ndarray:
    """Encode equal-length sequences to int8 (invalid/ambiguous = -1)."""
    n, L = len(seqs), len(seqs[0]) if seqs else 0
    arr = np.full((n, L), -1, dtype=np.int8)
    for i, s in enumerate(seqs):
        for j, ch in enumerate(s):
            v = BASES.get(ch)
            if v is not None:
                arr[i, j] = v
    return arr


def _best_scores(encoded: np.ndarray, pwm: np.ndarray, varidx: np.ndarray) -> np.ndarray:
    """Best log-odds score among sites containing varidx (per row), one strand.

    encoded: (N, L) int8 with values 0..3 (-1 tolerated -> zero contribution)
    pwm: (4, w) in ACGT-row order (as built by ReferenceBundle.motifs);
    varidx: (N,) position of the variant inside the sequence.
    """
    N, L = encoded.shape
    w = pwm.shape[1]
    best = np.full(N, -np.inf)
    for s in range(0, L - w + 1):
        # rows whose variant position lies within [s, s+w)
        hit = (varidx >= s) & (varidx < s + w)
        if not hit.any():
            continue
        sub = encoded[hit, s:s + w]
        # ambiguous/absent bases contribute 0 (masked one-hot)
        oh = np.zeros((sub.shape[0], w, 4), dtype=np.float64)
        m = sub >= 0
        rows = np.arange(sub.shape[0])[:, None]
        # oh[n, j, b] = 1 when sequence n has base b at position j of the site
        oh[rows, np.arange(w)[None, :], np.clip(sub, 0, 3)] = m.astype(float)
        # oh is (n, position, base); pwm is (base, position)
        sc = np.einsum("npb,bp->n", oh, pwm)
        idxs = np.nonzero(hit)[0]
        best[idxs] = np.maximum(best[idxs], sc)
    return best


def step_motif_effect(df: pd.DataFrame, *, bundle: ReferenceBundle,
                      config: AppConfig,
                      out_dir: Path) -> tuple[pd.DataFrame, FunnelStep, list[str]]:
    """Scan ref/alt windows for motif disruption; never drops rows."""
    rows_in = len(df)
    warnings: list[str] = []
    motifs = bundle.motifs()
    W = int(config.motif.window)
    min_strong = config.motif.min_strong
    min_drop = config.motif.min_drop

    df = df.copy()
    n_ref_mismatch = 0
    n_no_fasta = 0
    n_scanned = 0

    # precompute per-motif covering-site offset lists for ref and alt lengths
    motif_items = [(tf, m["pwm"], m["max_score"]) for tf, m in motifs.items()]

    out_cols = {k: np.zeros(len(df), dtype=np.float64) for k in
                ("motif_ref_rel", "motif_alt_rel", "motif_delta", "motif_disruption_score")}
    top_tf = np.full(len(df), None, dtype=object)
    effect = np.full(len(df), "none", dtype=object)

    for chrom, gi in df.groupby("chrom", sort=False).indices.items():
        if not bundle.has_fasta(chrom):
            n_no_fasta += len(gi)
            continue
        fa = bundle.fasta(chrom)
        clen = fa.get_reference_length(chrom)
        sub = df.iloc[gi]
        refs = sub["ref"].str.upper().to_numpy()
        alts = sub["alt"].str.upper().to_numpy()
        pos0 = sub["pos"].to_numpy(np.int64) - 1

        ref_seqs: list[str] = []
        alt_seqs: list[str] = []
        keep_local: list[int] = []
        for k, (p0, r, a) in enumerate(zip(pos0, refs, alts)):
            f_start = max(0, p0 - W)
            f_end = min(clen, p0 + len(r) + W)
            left_pad = (p0 - W) - f_start
            raw = fa.fetch(chrom, int(f_start), int(f_end)).upper()
            right_pad = (p0 + len(r) + W) - f_end
            window = ("N" * left_pad) + raw + ("N" * right_pad)
            if len(window) != 2 * W + len(r):
                continue
            if window[W:W + len(r)] != r:
                n_ref_mismatch += 1
                continue
            alt_window = window[:W] + a + window[W + len(r):]
            ref_seqs.append(window)
            alt_seqs.append(alt_window)
            keep_local.append(k)
        if not keep_local:
            continue
        n_scanned += len(keep_local)
        idx_local = np.asarray(keep_local)
        rows = gi[idx_local]

        n = len(rows)
        # the variant always starts at index W in both windows (SNV/indel anchor)
        varidx = np.full(n, W, dtype=np.int64)

        # window lengths vary with indels (2W + len(ref/alt)) -> group rows so
        # _encode sees uniform-length sequences; encodings are reused per motif
        ref_groups: list[tuple[np.ndarray, np.ndarray]] = []
        for L in sorted({len(s) for s in ref_seqs}):
            sel = np.array([k for k, s in enumerate(ref_seqs) if len(s) == L])
            ref_groups.append((sel, _encode([ref_seqs[k] for k in sel])))
        alt_groups: list[tuple[np.ndarray, np.ndarray]] = []
        for L in sorted({len(s) for s in alt_seqs}):
            sel = np.array([k for k, s in enumerate(alt_seqs) if len(s) == L])
            alt_groups.append((sel, _encode([alt_seqs[k] for k in sel])))

        # results for this chrom chunk
        best_delta = np.zeros(n)
        chosen_tf = np.full(n, None, dtype=object)
        chosen_ref = np.zeros(n)
        chosen_alt = np.zeros(n)
        chosen_eff = np.full(n, "none", dtype=object)
        any_effect = np.zeros(n, dtype=bool)

        for tf, pwm, max_score in motif_items:
            # rc-PWM: rows reversed (base complement) + columns reversed
            # (position flip); scored against the FORWARD sequence it detects
            # revcomp sites at the same variant-overlapping offsets
            pwm_rc = pwm[::-1, ::-1]
            # forward strand: forward sequence with the PWM
            r_f = np.empty(n)
            r_rc = np.empty(n)
            for sel, enc in ref_groups:
                r_f[sel] = _best_scores(enc, pwm, varidx[sel])
                r_rc[sel] = _best_scores(enc, pwm_rc, varidx[sel])
            a_f = np.empty(n)
            a_rc = np.empty(n)
            for sel, enc in alt_groups:
                a_f[sel] = _best_scores(enc, pwm, varidx[sel])
                a_rc[sel] = _best_scores(enc, pwm_rc, varidx[sel])
            r_best = np.maximum(r_f, r_rc) / max_score
            a_best = np.maximum(a_f, a_rc) / max_score
            r_rel = np.clip(r_best, 0, 1)
            a_rel = np.clip(a_best, 0, 1)
            delta = a_rel - r_rel

            is_broken = (r_rel >= min_strong) & ((r_rel - a_rel) >= min_drop)
            is_created = (a_rel >= min_strong) & ((a_rel - r_rel) >= min_drop)
            hit_tf = is_broken | is_created
            # stronger motif-disruption evidence wins; otherwise track max |delta|
            better = np.abs(delta) > np.abs(best_delta)
            take = better & (hit_tf | ~any_effect)
            # prioritize effectful variants over effectless ones
            upgrade = hit_tf & ~any_effect
            take = (upgrade) | (better & ~any_effect & ~hit_tf) | \
                   (better & any_effect & hit_tf)
            chosen_tf = np.where(take, tf, chosen_tf)
            chosen_ref = np.where(take, r_rel, chosen_ref)
            chosen_alt = np.where(take, a_rel, chosen_alt)
            chosen_eff = np.where(take,
                                  np.where(is_broken, "broken",
                                           np.where(is_created, "created", "none")),
                                  chosen_eff)
            best_delta = np.where(take, delta, best_delta)
            any_effect |= hit_tf

        out_cols["motif_ref_rel"][rows] = chosen_ref
        out_cols["motif_alt_rel"][rows] = chosen_alt
        out_cols["motif_delta"][rows] = best_delta
        out_cols["motif_disruption_score"][rows] = np.where(
            chosen_eff != "none", np.abs(best_delta), 0.0)
        top_tf[rows] = chosen_tf
        effect[rows] = chosen_eff

    for k, v in out_cols.items():
        df[k] = v
    df["top_motif_tf"] = top_tf
    df["motif_effect"] = effect
    df["motif_scanned"] = True
    df.loc[df["top_motif_tf"].isna(), "motif_scanned"] = False

    if n_no_fasta:
        raise FileNotFoundError(
            f"FASTA missing for {n_no_fasta} variants (chromosomes absent); run `make data`")
    if n_ref_mismatch:
        warnings.append(f"{n_ref_mismatch} variants skipped in s3 due to REF/FASTA mismatch "
                        "(possible build mismatch)")
    df.to_parquet(out_dir / "s3_motif.parquet", index=False)
    n_broken = int((df["motif_effect"] == "broken").sum())
    n_created = int((df["motif_effect"] == "created").sum())
    step = FunnelStep(
        step="s3", label="Motif disruption scan", rows_in=rows_in, rows_out=len(df),
        removed=0,
        note=f"scanned {n_scanned}; broken={n_broken}, created={n_created}; "
             f"{len(motifs)} PWMs; REF mismatches={n_ref_mismatch}")
    log.info("s3: scanned %d/%d (broken=%d created=%d mismatch=%d)",
             n_scanned, rows_in, n_broken, n_created, n_ref_mismatch)
    return df, step, warnings
