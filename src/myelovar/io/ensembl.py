"""Ensembl REST helpers: region variation overlap, LD proxies, population AF.

All requests are disk-cached (io.http) and rate-limited. Ensembl's main REST
server is GRCh38, so no liftover is needed for Mode B.
"""
from __future__ import annotations

import logging
from pathlib import Path

from myelovar.io.http import get_json

log = logging.getLogger(__name__)

REST = "https://rest.ensembl.org"
CACHE = Path("data/cache/http/ensembl")


def region_variations(chrom: str, start: int, end: int) -> list[dict]:
    """All dbSNP variations overlapping a GRCh38 region."""
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    url = f"{REST}/overlap/region/human/{chrom}:{start}-{end}?feature=variation"
    data = get_json(url, cache_dir=CACHE, timeout=90)
    if not isinstance(data, list):
        raise RuntimeError(f"unexpected Ensembl overlap payload: {str(data)[:200]}")
    return data


def ld_proxies(rsid: str, population: str, r2: float) -> list[dict]:
    """Variants in LD (r2 >= threshold) with rsid in a 1000G population."""
    url = f"{REST}/ld/human/{rsid}/{population}?r2={r2}"
    try:
        data = get_json(url, cache_dir=CACHE, timeout=120)
    except RuntimeError as exc:
        log.warning("LD query failed for %s: %s", rsid, exc)
        return []
    return data if isinstance(data, list) else []


def variation_detail(rsid: str, with_pops: bool = False) -> dict:
    """Variation record incl. GRCh38 mapping (and population AFs when asked)."""
    suffix = "?pops=1" if with_pops else ""
    url = f"{REST}/variation/human/{rsid}{suffix}"
    return get_json(url, cache_dir=CACHE, timeout=60)


def eur_af(detail: dict, population: str = "1000GENOMES:phase_3:EUR",
           alt: str | None = None) -> float | None:
    """Extract a population allele frequency (for the given ALT when provided)."""
    pops = detail.get("populations") or []
    vals = []
    for p in pops:
        if p.get("population") == population and p.get("allele") not in (None, "C", "G"):
            pass
    for p in pops:
        if p.get("population") != population:
            continue
        allele = (p.get("allele") or "").upper()
        freq = p.get("frequency")
        if freq is None:
            continue
        if alt is None:
            if allele not in ("A", "C", "G", "T"):  # major/other buckets excluded when possible
                continue
            vals.append(float(freq))
        elif allele == alt.upper():
            return float(freq)
    return max(vals) if vals else None
