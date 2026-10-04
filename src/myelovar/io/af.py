"""Bulk allele-frequency lookup from the UCSC dbSNP-common track (hg38).

Why this and not gnomAD: the gnomAD sites VCFs are 35-53 GB *per chromosome*
(verified at runtime), infeasible on this host. The UCSC snp151Common table
carries real submitted allele frequencies (`alleleFreqs`) for common variants
(MAF >= 1% in a surveyed population), which implements the spec's
"drop AF > 0.01" filter with real public data. Shortlist variants additionally
get gnomAD AF from Ensembl VEP / the gnomAD browser API (see io/ensembl.py).

Matching: (chrom, pos, ref, alt) against dbSNP refUCSC/refNCBI + observed
alleles. If a variant matches the common track but no numeric frequency can be
parsed for its ALT, it is still treated as common (track membership is itself a
common-variant criterion) with reason `in_dbSNP_common_track`.
"""
from __future__ import annotations

import gzip
import hashlib
import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# snp151Common column indices (from the published CREATE TABLE, verified at runtime)
C_CHROM, C_START, C_NAME, C_REFNCBI, C_REFUCSC = 1, 2, 4, 7, 8
C_OBSERVED, C_AVHET, C_FREQCOUNT, C_ALLELES, C_FREQS = 9, 13, 21, 22, 24


def _split_alleles(observed: str, alleles: str) -> list[str]:
    out: list[str] = []
    for blob in (observed, alleles):
        if not blob or blob == "?":
            continue
        for part in blob.replace("/", ",").replace(";", ",").split(","):
            part = part.strip().upper()
            if part and part != "?" and part != "-":
                out.append(part)
    return out


def _alt_freq(fields: list[str], alt: str, ref: str) -> tuple[float | None, bool]:
    """Return (freq for ALT, parse_ok) from a dbSNP-common row."""
    observed = fields[C_OBSERVED] if len(fields) > C_OBSERVED else ""
    alleles = fields[C_ALLELES] if len(fields) > C_ALLELES else ""
    freqs_raw = fields[C_FREQS] if len(fields) > C_FREQS else ""
    count = fields[C_FREQCOUNT] if len(fields) > C_FREQCOUNT else "0"
    alt_u = alt.upper()
    try:
        n = int(count)
    except ValueError:
        n = 0
    # preferred: allele list paired with alleleFreqs
    if n > 0 and freqs_raw:
        allele_list = [a.strip().upper() for a in alleles.split(",") if a.strip()] \
            if alleles and alleles != "?" else []
        freq_list: list[str] = [f.strip() for f in freqs_raw.split(",") if f.strip()]
        if len(allele_list) == len(freq_list) == n:
            for a, f in zip(allele_list, freq_list):
                if a == alt_u:
                    try:
                        return float(f), True
                    except ValueError:
                        break
        # fallback pairing: observed alleles with freqs
        obs_list = [a.strip().upper() for a in observed.replace("/", ",").split(",") if a.strip()]
        if len(obs_list) == len(freq_list) == n:
            for a, f in zip(obs_list, freq_list):
                if a == alt_u:
                    try:
                        return float(f), True
                    except ValueError:
                        break
    return None, False


def extract_af(variants: pd.DataFrame, dbsnp_path: Path | str,
               cache_dir: Path | str, cache_key: str) -> pd.DataFrame:
    """Return per-variant (af, rsid, in_common_track) columns for `variants`.

    Results are cached to ``cache_dir/af/<cache_key>.parquet``.
    """
    cache = Path(cache_dir) / "af"
    cache.mkdir(parents=True, exist_ok=True)
    cache_file = cache / f"{cache_key}.parquet"
    if cache_file.exists():
        cached = pd.read_parquet(cache_file)
        if len(cached) == len(variants) and (cached["variant_id"].values == variants["variant_id"].values).all():
            log.info("AF cache hit: %s (%d rows)", cache_file, len(cached))
            return cached
        log.warning("AF cache %s mismatched input; recomputing", cache_file)

    dbsnp_path = Path(dbsnp_path)
    if not dbsnp_path.exists():
        raise FileNotFoundError(
            f"dbSNP common track not found at {dbsnp_path}; run `make data` first")

    # key positions for an O(1) probe
    wanted: dict[tuple[str, int], list[tuple[str, str, str]]] = {}
    for vid, chrom, pos, ref, alt in zip(variants["variant_id"], variants["chrom"],
                                         variants["pos"], variants["ref"], variants["alt"]):
        wanted.setdefault((chrom, int(pos)), []).append((ref.upper(), alt.upper(), vid))

    found: dict[str, tuple[float | None, str, bool]] = {}
    n_rows = 0
    with gzip.open(dbsnp_path, "rt") as fh:
        for line in fh:
            n_rows += 1
            # cheap pre-filter before splitting
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            key = (f[C_CHROM], int(f[C_START]) + 1)
            cands = wanted.get(key)
            if not cands:
                continue
            ref_ucsc = (f[C_REFUCSC] or "").upper()
            ref_ncbi = (f[C_REFNCBI] or "").upper()
            pool = set(_split_alleles(f[C_OBSERVED], f[C_ALLELES] if len(f) > C_ALLELES else ""))
            for ref, alt, vid in cands:
                if ref not in (ref_ucsc, ref_ncbi):
                    # indels often differ in representation; require alt presence anyway
                    if alt not in pool:
                        continue
                elif alt not in pool and len(ref) == 1:
                    # SNV whose observed list should contain the alt; skip if absent
                    continue
                freq, parsed = _alt_freq(f, alt, ref)
                found[vid] = (freq, f[C_NAME], True)
    n_rows_scanned = n_rows

    out = pd.DataFrame({
        "variant_id": variants["variant_id"].to_numpy(),
        "af": [found.get(v, (None, None, False))[0] for v in variants["variant_id"]],
        "rsid": [found.get(v, (None, None, False))[1] for v in variants["variant_id"]],
        "in_common_track": [found.get(v, (None, None, False))[2] for v in variants["variant_id"]],
    })
    out["af"] = pd.to_numeric(out["af"], errors="coerce")
    log.info("AF extraction: scanned %d dbSNP-common rows; matched %d/%d variants",
             n_rows_scanned, len(found), len(variants))
    out.to_parquet(cache_file, index=False)
    return out
