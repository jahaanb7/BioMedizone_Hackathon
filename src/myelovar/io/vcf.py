"""VCF streaming loader (pysam/htslib) -> tidy DataFrame (s0).

Normalization rules:
  * chromosome names gain a 'chr' prefix if missing;
  * multiallelic records are split into one row per ALT allele;
  * variant_id = chr:pos:ref:alt (post-split);
  * DP/GQ are read from FORMAT when present, else kept None (logged by s1).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import pandas as pd
import pysam

log = logging.getLogger(__name__)

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def normalize_chrom(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def load_vcf(path: Path | str, *, autosomes_only: bool = True,
             max_records: Optional[int] = None,
             chrom_subset: Optional[list[str]] = None) -> pd.DataFrame:
    """Stream a VCF/VCF.gz into a DataFrame; one row per (site, ALT allele)."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"input VCF not found: {path}")
    if path.suffix not in (".vcf", ".gz", ".bgz") and not path.name.endswith(".vcf.gz"):
        raise ValueError(f"input must be .vcf or .vcf.gz, got: {path.name}")

    rows: list[dict] = []
    n_records = 0
    n_split = 0
    has_dp = has_gq = False
    allowed = set(chrom_subset) if chrom_subset else (
        set(AUTOSOMES) if autosomes_only else None)

    vf = pysam.VariantFile(str(path))
    sample = next(iter(vf.header.samples), None)
    for rec in vf:
        n_records += 1
        chrom = normalize_chrom(rec.chrom)
        if allowed is not None and chrom not in allowed:
            continue
        filt = list(rec.filter.keys()) if rec.filter else []
        filt_str = ";".join(filt) if filt else "."
        fmt = rec.samples[sample] if sample else None
        dp = gq = None
        if fmt is not None:
            if "DP" in fmt:
                dp = fmt.get("DP")
                has_dp = has_dp or dp is not None
            if "GQ" in fmt:
                gq = fmt.get("GQ")
                has_gq = has_gq or gq is not None
            gt = fmt.get("GT")
            gt_str = "/".join("." if a is None else str(a) for a in gt) if gt else "."
        else:
            gt_str = "."
        alts = rec.alts or ()
        if len(alts) > 1:
            n_split += 1
        for alt in alts:
            rows.append({
                "chrom": chrom, "pos": rec.pos, "ref": rec.ref, "alt": alt,
                "variant_id": f"{chrom}:{rec.pos}:{rec.ref}:{alt}",
                "id": rec.id or ".",
                "qual": rec.qual, "filter": filt_str,
                "dp": int(dp) if dp is not None else None,
                "gq": int(gq) if gq is not None else None,
                "gt": gt_str,
            })
        if max_records is not None and n_records >= max_records:
            break
    vf.close()

    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError(f"no variants loaded from {path} (check chromosome filter)")
    df = df.drop_duplicates(subset="variant_id", keep="first").reset_index(drop=True)
    log.info("s0 loaded %d records -> %d variant rows (multiallelic splits=%d; "
             "DP present=%s, GQ present=%s)", n_records, len(df), n_split, has_dp, has_gq)
    return df
