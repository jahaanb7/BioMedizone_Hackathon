"""ReferenceBundle: one lazily-loaded, cached handle to all reference data.

Heavy data (track interval indexes, GENCODE tables, motif matrices, DepMap
summary, ABC predictions) loads once per process and is precomputed on disk by
scripts/build_cache.py, so the API server starts fast.
"""
from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from myelovar.config import AppConfig
from myelovar.io.af import extract_af
from myelovar.io.bed import IntervalIndex, load_bed, open_text
from myelovar.io.gtf import parse_gencode

log = logging.getLogger(__name__)


class ReferenceBundle:
    """Lazy loader for tracks, annotation, motifs, DepMap, ABC and loci."""

    def __init__(self, config: AppConfig, root: Path | str = "."):
        self.config = config
        self.root = Path(root)
        self.data = self.root / config.paths.data_dir
        self.cache = self.root / config.paths.cache_dir
        self._manifest: Optional[dict] = None
        self._registry: Optional[dict] = None
        self._track_idx: dict[str, IntervalIndex] = {}
        self._gencode: Optional[dict] = None
        self._cds_idx: Optional[IntervalIndex] = None
        self._exon_idx: Optional[IntervalIndex] = None
        self._motifs: Optional[dict] = None
        self._depmap: Optional[pd.DataFrame] = None
        self._mm_genes: Optional[pd.DataFrame] = None
        self._abc: Optional[tuple[pd.DataFrame, dict]] = None
        self._loci: Optional[list[dict]] = None
        self._fastas: dict[str, Any] = {}
        self._blacklist: Optional[IntervalIndex] = None
        self._ccre: Optional[IntervalIndex] = None

    # ---------------------------------------------------------- manifest ---
    def manifest(self) -> dict:
        if self._manifest is None:
            path = self.root / self.config.paths.manifest
            if path.exists():
                self._manifest = json.loads(path.read_text())
            else:
                log.warning("data manifest missing at %s (run `make data`)", path)
                self._manifest = {"generated_at": "", "manifest_hash": "",
                                  "genome_build": "GRCh38", "resources": []}
        return self._manifest

    @property
    def manifest_hash(self) -> str:
        return str(self.manifest().get("manifest_hash", ""))

    def manifest_entry(self, resource_id: str) -> Optional[dict]:
        for r in self.manifest().get("resources", []):
            if r.get("id") == resource_id:
                return r
        return None

    # ------------------------------------------------------------- tracks ---
    def track_registry(self) -> dict:
        if self._registry is None:
            path = self.root / self.config.paths.track_registry
            if not path.exists():
                raise FileNotFoundError(
                    f"track registry {path} missing; run `make data` (and `scripts/build_cache.py`)")
            self._registry = json.loads(path.read_text())
        return self._registry

    def peak_tracks(self) -> list[dict]:
        return [t for t in self.track_registry()["tracks"] if t.get("kind") == "peaks"]

    def signal_tracks(self) -> list[dict]:
        return [t for t in self.track_registry()["tracks"] if t.get("kind") == "signal"]

    def track_index(self, name: str) -> IntervalIndex:
        if name in self._track_idx:
            return self._track_idx[name]
        entry = next((t for t in self.track_registry()["tracks"] if t["name"] == name), None)
        if entry is None:
            raise KeyError(f"unknown track: {name}")
        cache = self.cache / "tracks" / f"{name}.parquet"
        if cache.exists():
            df = pd.read_parquet(cache)
        else:
            df, _ = load_bed(self.root / entry["path"])
            cache.parent.mkdir(parents=True, exist_ok=True)
            df.to_parquet(cache, index=False)
        idx = IntervalIndex(df, merge=True)
        self._track_idx[name] = idx
        log.info("track index '%s': %d merged intervals", name, len(idx))
        return idx

    # ------------------------------------------------------- annotation ---
    def gencode(self) -> dict[str, pd.DataFrame]:
        if self._gencode is None:
            gtfs = sorted((self.data / "reference").glob("gencode.v*.basic.annotation.gtf.gz"))
            if not gtfs:
                raise FileNotFoundError("GENCODE GTF not found in data/reference; run `make data`")
            self._gencode = parse_gencode(gtfs[-1], self.cache)
        return self._gencode

    def cds_index(self) -> IntervalIndex:
        if self._cds_idx is None:
            self._cds_idx = IntervalIndex(self.gencode()["cds"], merge=True)
        return self._cds_idx

    def exon_index(self) -> IntervalIndex:
        if self._exon_idx is None:
            self._exon_idx = IntervalIndex(self.gencode()["exons"], merge=True)
        return self._exon_idx

    def genes(self) -> pd.DataFrame:
        return self.gencode()["genes"]

    def nearest_gene(self, chrom: str, pos0: int, max_scan: int = 400) -> tuple[Optional[str], Optional[int]]:
        """(gene, distance_bp) for the gene span nearest pos0 (0 if inside)."""
        g = self.genes()
        sub = g[g["chrom"] == chrom]
        if sub.empty:
            return None, None
        starts = sub["start"].to_numpy()
        i = int(np.searchsorted(starts, pos0, side="right")) - 1
        best_gene, best_dist = None, None
        for j in range(max(0, i - max_scan), min(len(sub), i + 2)):
            if j < 0:
                continue
            row = sub.iloc[j]
            s, e = int(row["start"]), int(row["end"])
            d = 0 if s <= pos0 < e else (pos0 - e if pos0 >= e else s - pos0)
            if best_dist is None or d < best_dist:
                best_dist, best_gene = d, str(row["gene"])
        return best_gene, best_dist

    def genes_tss_within(self, chrom: str, pos0: int, window: int) -> pd.DataFrame:
        """Genes whose TSS lies within +/- window bp of pos0."""
        g = self.genes()
        sub = g[g["chrom"] == chrom]
        if sub.empty:
            return sub
        tss = sub["tss"].to_numpy()
        lo = int(np.searchsorted(tss, pos0 - window, side="left"))
        hi = int(np.searchsorted(tss, pos0 + window, side="right"))
        return sub.iloc[lo:hi]

    # --------------------------------------------------------- blacklist ---
    def blacklist(self) -> IntervalIndex:
        if self._blacklist is None:
            path = self.data / "reference/hg38-blacklist.v2.bed.gz"
            if not path.exists():
                raise FileNotFoundError(f"blacklist {path} missing; run `make data`")
            df, _ = load_bed(path)
            self._blacklist = IntervalIndex(df, merge=True)
        return self._blacklist

    def ccre(self) -> IntervalIndex:
        if self._ccre is None:
            path = self.data / "tracks/ccre/ccre_v4_grch38.bed"
            if not path.exists():
                raise FileNotFoundError(f"cCRE file {path} missing; run `make data`")
            df, _ = load_bed(path)
            self._ccre = IntervalIndex(df, merge=True)
        return self._ccre

    # ------------------------------------------------------------ motifs ---
    def motifs(self) -> dict[str, dict]:
        """JASPAR PWMs as log-odds matrices (dict keyed by TF symbol)."""
        if self._motifs is None:
            path = self.cache / "jaspar/matrices.json"
            if not path.exists():
                raise FileNotFoundError(f"JASPAR matrices {path} missing; run `make data`")
            raw = json.loads(path.read_text())["matrices"]
            motifs: dict[str, dict] = {}
            for tf, m in raw.items():
                pfm = m["pfm"]
                counts = np.vstack([pfm[b] for b in "ACGT"]).astype(np.float64)
                total = counts.sum(axis=0, keepdims=True)
                total[total == 0] = 1.0
                probs = (counts + 0.25) / (total + 1.0)   # Laplace-ish pseudocount
                bg = 0.25
                pwm = np.log2(probs / bg)
                motifs[tf] = {"matrix_id": m["matrix_id"], "pwm": pwm,
                              "max_score": float(pwm.max(axis=0).sum()),
                              "min_score": float(pwm.min(axis=0).sum())}
            self._motifs = motifs
            log.info("loaded %d motif matrices", len(motifs))
        return self._motifs

    # ----------------------------------------------------------- depmap ---
    def depmap(self) -> pd.DataFrame:
        if self._depmap is None:
            path = self.cache / "depmap/depmap_summary.parquet"
            if not path.exists():
                raise FileNotFoundError(
                    f"DepMap summary {path} missing; run `scripts/build_cache.py`")
            self._depmap = pd.read_parquet(path)
        return self._depmap

    def mm_genes(self) -> pd.DataFrame:
        if self._mm_genes is None:
            path = self.data / "reference/mm_genes.tsv"
            if path.exists():
                self._mm_genes = pd.read_csv(path, sep="\t")
            else:
                log.warning("mm_genes.tsv missing (run `make data`); known-gene flag disabled")
                self._mm_genes = pd.DataFrame(columns=["gene", "source"])
        return self._mm_genes

    # -------------------------------------------------------------- abc ---
    def abc(self) -> tuple[pd.DataFrame, dict]:
        if self._abc is None:
            path = self.cache / "abc/abc_predictions.parquet"
            meta_path = self.cache / "abc/abc_meta.json"
            if not path.exists():
                raise FileNotFoundError(f"ABC cache {path} missing; run `scripts/build_cache.py`")
            meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
            self._abc = (pd.read_parquet(path), meta)
        return self._abc

    # ------------------------------------------------------------- loci ---
    def gwas_loci(self) -> list[dict]:
        """Myeloma GWAS loci (lead SNPs +/- window) for in_mm_gwas_locus."""
        if self._loci is None:
            path = self.cache / "gwas/loci.json"
            if path.exists():
                self._loci = json.loads(path.read_text())
            else:
                log.warning("gwas/loci.json missing; in_mm_gwas_locus will be False "
                            "(run `myelovar run --mode B` or scripts/build_cache.py)")
                self._loci = []
        return self._loci

    # ------------------------------------------------------------- fasta ---
    def fasta(self, chrom: str):
        if chrom not in self._fastas:
            path = self.data / f"reference/genome/{chrom}.fa.bgz"
            if not path.exists():
                raise FileNotFoundError(
                    f"reference FASTA for {chrom} missing ({path}); run `make data`")
            import pysam
            self._fastas[chrom] = pysam.FastaFile(str(path))
        return self._fastas[chrom]

    def has_fasta(self, chrom: str) -> bool:
        return (self.data / f"reference/genome/{chrom}.fa.bgz").exists()

    # ---------------------------------------------------------------- AF ---
    def af(self, variants: pd.DataFrame) -> pd.DataFrame:
        """dbSNP-common AF lookup, cached per input variant set."""
        key = hashlib.sha1(",".join(variants["variant_id"]).encode()).hexdigest()[:20]
        dbsnp = self.data / "raw/dbsnp/snp151Common.txt.gz"
        return extract_af(variants, dbsnp, self.cache, key)
