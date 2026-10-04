"""GENCODE GTF parsing: CDS exons, gene spans, TSS table (cached as parquet)."""
from __future__ import annotations

import gzip
import logging
import re
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

_ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


def parse_gencode(gtf_path: Path | str, cache_dir: Path | str) -> dict[str, pd.DataFrame]:
    """Parse CDS/exon/gene records; returns {'cds','exons','genes'} frames.

    Cached per-GTF-filename under cache_dir/gencode/.
    """
    gtf_path = Path(gtf_path)
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF not found: {gtf_path}; run `make data`")
    cache = Path(cache_dir) / "gencode" / gtf_path.name.replace(".gz", "")
    cache.mkdir(parents=True, exist_ok=True)
    paths = {k: cache / f"{k}.parquet" for k in ("cds", "exons", "genes")}
    if all(p.exists() for p in paths.values()):
        return {k: pd.read_parquet(p) for k, p in paths.items()}

    cds_rows, exon_rows, gene_rows = [], [], []
    n = 0
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feature = f[2]
            if feature not in ("CDS", "exon", "gene"):
                continue
            chrom = f[0] if f[0].startswith("chr") else f"chr{f[0]}"
            start, end = int(f[3]) - 1, int(f[4])  # GTF 1-based inclusive -> 0-based half-open
            strand = f[6]
            attrs = dict(_ATTR_RE.findall(f[8]))
            if feature == "gene":
                gene_rows.append((chrom, start, end, attrs.get("gene_name", attrs.get("gene_id", "?")),
                                  attrs.get("gene_type", attrs.get("gene_biotype", "?")), strand))
            elif feature == "exon":
                exon_rows.append((chrom, start, end))
            else:
                cds_rows.append((chrom, start, end))
            n += 1

    genes = pd.DataFrame(gene_rows, columns=["chrom", "start", "end", "gene", "gene_type", "strand"])
    genes["tss"] = genes.apply(lambda r: r["start"] if r["strand"] == "+" else r["end"] - 1, axis=1)
    cds = pd.DataFrame(cds_rows, columns=["chrom", "start", "end"])
    exons = pd.DataFrame(exon_rows, columns=["chrom", "start", "end"])
    if genes.empty or cds.empty:
        raise RuntimeError(f"GENCODE parse produced empty tables from {gtf_path}")
    genes = genes.drop_duplicates(subset=["chrom", "gene", "tss"])
    genes.to_parquet(paths["genes"], index=False)
    cds.to_parquet(paths["cds"], index=False)
    exons.to_parquet(paths["exons"], index=False)
    log.info("GENCODE parsed: %d genes, %d CDS, %d exons (from %d feature lines)",
             len(genes), len(cds), len(exons), n)
    return {"cds": cds, "exons": exons, "genes": genes}
