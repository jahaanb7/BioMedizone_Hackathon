"""ABC enhancer-gene predictions (Nasser et al. Nature 2021): cache a
plasma/B-lineage biosample subset from the all-biosamples file.

The chosen CellType and whether it is a proxy (not an actual plasma-cell
biosample) are recorded in the cache metadata and surfaced in the UI/manifest.
"""
from __future__ import annotations

import gzip
import json
import logging
import re
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

# priority order; first class with a hit wins, proxy flag set when != plasma
_PREFS = [
    (r"plasma", False),
    (r"bone.?marrow", True),
    (r"\bB[-_ ]?cell|Bcell|Bcells|naive.?B|memory.?B|GC.?B|tonsil|lymphoblastoid|GM12878", True),
    (r"leukemia|myeloma|lymphoma", True),
]


def _choose_celltypes(counts: dict[str, int]) -> tuple[list[str], bool]:
    for pat, proxy in _PREFS:
        hits = sorted([ct for ct in counts if re.search(pat, ct, re.I)])
        if hits:
            return hits, (proxy or not re.search(r"plasma", pat, re.I))
    # last resort: largest biosample
    top = max(counts, key=counts.get)
    return [top], True


def build_abc_cache(raw_path: Path | str, cache_dir: Path | str, *,
                    force: bool = False) -> dict:
    """Filter the ABC all-predictions file to the best plasma/B-lineage CellType(s).

    Returns metadata dict; writes <cache_dir>/abc/abc_predictions.parquet.
    """
    raw_path = Path(raw_path)
    if not raw_path.exists():
        raise FileNotFoundError(f"ABC predictions not found at {raw_path}; run `make data`")
    cache = Path(cache_dir) / "abc"
    cache.mkdir(parents=True, exist_ok=True)
    parquet = cache / "abc_predictions.parquet"
    meta_path = cache / "abc_meta.json"
    if parquet.exists() and meta_path.exists() and not force:
        return json.loads(meta_path.read_text())

    # pass 1: cell type inventory
    counts: dict[str, int] = {}
    with gzip.open(raw_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        need = ["chr", "start", "end", "TargetGene", "TargetGeneTSS", "ABC.Score", "CellType", "distance"]
        missing = [c for c in need if c not in idx]
        if missing:
            raise RuntimeError(f"ABC header missing {missing}: {header}")
        for line in fh:
            ct = line.rstrip("\n").rsplit("\t", 1)[-1]
            counts[ct] = counts.get(ct, 0) + 1
    chosen, is_proxy = _choose_celltypes(counts)
    chosen_set = set(chosen)
    log.info("ABC: %d biosamples; chosen %s (proxy=%s)", len(counts), chosen, is_proxy)

    # pass 2: extract rows for the chosen cell types
    rows = []
    with gzip.open(raw_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            # cheap guard before full split
            if not any(ct in line for ct in chosen_set):
                continue
            f = line.rstrip("\n").split("\t")
            if f[idx["CellType"]] not in chosen_set:
                continue
            rows.append((f[idx["chr"]], int(f[idx["start"]]), int(f[idx["end"]]),
                         f[idx["TargetGene"]], int(float(f[idx["TargetGeneTSS"]])),
                         float(f[idx["ABC.Score"]]), f[idx["CellType"]],
                         int(float(f[idx["distance"]]))))
    if not rows:
        raise RuntimeError(f"no ABC rows for chosen cell types {chosen}")
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene", "tss",
                                     "abc_score", "celltype", "distance"])
    df.to_parquet(parquet, index=False)
    meta = {"cell_types": chosen, "proxy": is_proxy,
            "n_biosamples_total": len(counts), "n_rows": len(df),
            "source": str(raw_path),
            "note": ("chosen biosample is NOT a plasma-cell sample; used as a labelled proxy"
                     if is_proxy else "plasma-cell biosample")}
    meta_path.write_text(json.dumps(meta, indent=2))
    log.info("ABC cache: %d predictions for %s", len(df), chosen)
    return meta
