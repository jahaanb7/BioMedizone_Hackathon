"""s0: load input variants (Mode A VCF or Mode B GWAS-locus table)."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from myelovar.config import AppConfig
from myelovar.io.vcf import load_vcf
from myelovar.schemas import FunnelStep

log = logging.getLogger(__name__)


def step_load_input(config: AppConfig, *, input_path: Optional[str] = None,
                    max_records: Optional[int] = None) -> tuple[pd.DataFrame, FunnelStep]:
    """Load variants for the configured mode.

    Mode A: stream a VCF/VCF.gz (autosomes by default; multiallelics split).
    Mode B: build the locus variant table from the cached GWAS Catalog +
            Ensembl data (io.gwas.build_mode_b_variants).
    """
    mode = config.mode.upper()
    if mode == "B":
        from myelovar.io.gwas import build_mode_b_variants
        df = build_mode_b_variants(config, force=False)
        desc = f"Mode B GWAS-locus table ({len(df)} candidate variants, {df['locus_id'].nunique()} loci)"
    elif mode == "A":
        path = input_path or config.input
        if not path:
            # default demo input: the GIAB benchmark genome (recorded in manifest)
            path = "data/raw/giab/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
        df = load_vcf(path, autosomes_only=config.filter.autosomes_only,
                      max_records=max_records)
        desc = f"Mode A VCF {Path(path).name} ({len(df)} variant rows)"
    elif mode == "C":
        from myelovar.io.gwas import build_mode_c_table  # stretch: DepMap somatic
        df = build_mode_c_table(config)
        desc = f"Mode C DepMap somatic ({len(df)} variants)"
    else:
        raise ValueError(f"unknown mode {config.mode!r}; expected A, B or C")

    step = FunnelStep(step="s0", label="Load input", rows_in=len(df),
                      rows_out=len(df), removed=0,
                      note=desc)
    log.info("s0: %s", desc)
    return df, step
