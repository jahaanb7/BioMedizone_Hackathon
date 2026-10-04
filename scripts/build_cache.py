#!/usr/bin/env python3
"""Build every precomputed cache the pipeline/API reads at request time.

Outputs:
  data/cache/gencode/*            parsed GENCODE tables (parquet)
  data/cache/af/<hash>.parquet    dbSNP-common AF for the default Mode A input
  data/cache/abc/*                ABC predictions filtered to chosen biosample
  data/cache/depmap/*             myeloma-vs-other gene effect summary
  data/cache/gwas/loci.json       merged myeloma GWAS loci
  data/cache/tracks/*.parquet     per-track interval indexes
  data/tracks/super_enhancer.bed  ROSE-style super-enhancer calls (MM.1S)
  (appends derived + GEO tracks to data/tracks/tracks.json)

Run: python scripts/build_cache.py [--skip-af]
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("build_cache")

ROOT = Path(__file__).resolve().parent.parent


def _load_registry(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text())
    raise FileNotFoundError(f"{path} missing - run `make data` first")


def _save_registry(path: Path, reg: dict) -> None:
    path.write_text(json.dumps(reg, indent=2))


def _upsert_track(reg: dict, entry: dict) -> None:
    tracks = reg.setdefault("tracks", [])
    tracks[:] = [t for t in tracks if t.get("name") != entry["name"]]
    tracks.append(entry)


# ------------------------------------------------------------ GENCODE ---
def build_gencode(cfg) -> None:
    from myelovar.io.gtf import parse_gencode
    gtfs = sorted((ROOT / "data/reference").glob("gencode.v*.basic.annotation.gtf.gz"))
    if not gtfs:
        raise FileNotFoundError("GENCODE GTF missing; run `make data`")
    parse_gencode(gtfs[-1], ROOT / cfg.paths.cache_dir)
    log.info("GENCODE cache built")


# --------------------------------------------------------------- ABC ---
def build_abc(cfg) -> None:
    from myelovar.io.abc import build_abc_cache
    raw = ROOT / "data/raw/abc/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
    meta = build_abc_cache(raw, ROOT / cfg.paths.cache_dir)
    log.info("ABC cache: %s", meta)


# ------------------------------------------------------------ DepMap ---
def build_depmap(cfg) -> None:
    from myelovar.io.depmap import build_depmap_summary
    summary = build_depmap_summary(
        ROOT / cfg.paths.cache_dir,
        mm_effect_threshold=cfg.depmap.mm_effect_threshold,
        selectivity_margin=cfg.depmap.selectivity_margin,
        expressed_pct_threshold=cfg.depmap.expressed_pct_threshold)
    log.info("DepMap summary: %d genes, %d selective", len(summary),
             int(summary["depmap_mm_selective"].sum()))


# -------------------------------------------------------------- GWAS ---
def build_gwas_loci(cfg) -> None:
    from myelovar.io.gwas import build_loci
    loci = build_loci(cfg, ROOT / "data/cache/gwas")
    log.info("GWAS loci: %d", len(loci))


# ------------------------------------------------------ super-enhancers ---
def _read_peaks_with_signal(path: Path) -> "tuple[np.ndarray, np.ndarray, np.ndarray]":
    """Read peak BED(narrowPeak) with per-interval signal (col 7 when present)."""
    import pandas as pd
    rows = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line or line[0] in "#t" or line.startswith("browser"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            sig = None
            if len(f) >= 7:
                try:
                    sig = float(f[6])
                except ValueError:
                    sig = None
            if sig is None and len(f) >= 5:
                try:
                    sig = float(f[4])
                except ValueError:
                    sig = None
            rows.append((f[0], int(f[1]), int(f[2]), sig if sig is not None else 1.0))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "signal"])
    return df


def build_super_enhancer(cfg, reg: dict) -> None:
    """ROSE-style: stitch H3K27ac peaks within 12.5 kb, exclude promoter-only
    elements (TSS +/-2 kb), rank by summed H3K27ac signal, cut at the ranked
    curve's inflection point (max distance from first->last chord)."""
    from myelovar.io.bed import IntervalIndex
    from myelovar.reference import ReferenceBundle

    peak_tracks = [t for t in reg.get("tracks", [])
                   if t.get("category") == "h3k27ac" and t.get("kind") == "peaks"]
    if not peak_tracks:
        log.warning("no H3K27ac peak track; skipping super-enhancer build")
        return
    primary = next((t for t in peak_tracks if t["experiment"].startswith("ENC")), peak_tracks[0])
    peaks = _read_peaks_with_signal(ROOT / primary["path"])
    log.info("ROSE: %d H3K27ac peaks from %s", len(peaks), primary["name"])

    stitched: list[tuple[str, int, int, float]] = []
    for chrom, sub in peaks.groupby("chrom", sort=False):
        sub = sub.sort_values("start")
        cs, ce, sig = int(sub["start"].iloc[0]), int(sub["end"].iloc[0]), float(sub["signal"].iloc[0])
        for s, e, sg in zip(sub["start"].iloc[1:], sub["end"].iloc[1:], sub["signal"].iloc[1:]):
            if int(s) - int(ce) <= 12_500:
                ce = max(ce, int(e))
                sig += float(sg)
            else:
                stitched.append((chrom, cs, ce, sig))
                cs, ce, sig = int(s), int(e), float(sg)
        stitched.append((chrom, cs, ce, sig))

    # promoter exclusion (elements overlapping any GENCODE TSS +/-2 kb)
    bundle = ReferenceBundle(cfg, root=ROOT)
    genes = bundle.genes()
    tss_by_chrom: dict[str, np.ndarray] = {
        c: np.sort(g["tss"].to_numpy(np.int64))
        for c, g in genes.groupby("chrom")}
    enhancer_elements = []
    n_promoter = 0
    for chrom, cs, ce, sig in stitched:
        arr = tss_by_chrom.get(chrom)
        if arr is not None and len(arr):
            lo = int(np.searchsorted(arr, cs - 2000, side="left"))
            hi = int(np.searchsorted(arr, ce + 2000, side="right"))
            if hi > lo:
                n_promoter += 1
                continue
        enhancer_elements.append((chrom, cs, ce, sig))
    log.info("ROSE: %d stitched, %d promoter-only excluded, %d ranked",
             len(stitched), n_promoter, len(enhancer_elements))

    # signal: prefer bigWig sum (pyBigWig), fall back to summed peak signals
    method = "summed narrowPeak signalValue"
    sig_tracks = [t for t in reg.get("tracks", [])
                  if t.get("category") == "h3k27ac" and t.get("kind") == "signal"
                  and t.get("experiment") == primary["experiment"]]
    if sig_tracks:
        try:
            import pyBigWig
            bw = pyBigWig.open(str(ROOT / sig_tracks[0]["path"]))
            if bw.isBigWig():
                for i, (chrom, cs, ce, _sig) in enumerate(enhancer_elements):
                    try:
                        v = bw.values(chrom, cs, ce, numpy=True)
                        total = float(np.nansum(v))
                        if np.isfinite(total) and total > 0:
                            enhancer_elements[i] = (chrom, cs, ce, total)
                    except Exception:  # noqa: BLE001
                        continue
                method = f"pyBigWig sum over {sig_tracks[0]['name']} ({sig_tracks[0].get('output_type')})"
            bw.close()
        except Exception as exc:  # noqa: BLE001
            log.warning("pyBigWig signal unavailable (%s); using peak signalValue", exc)

    if len(enhancer_elements) < 3:
        log.warning("too few enhancer elements for SE call (%d); skipping",
                    len(enhancer_elements))
        return
    ranked = sorted(enhancer_elements, key=lambda x: -x[3])
    y = np.array([r[3] for r in ranked])
    n = len(y)
    x = np.arange(n)
    # distance from chord (0,y0)-(n-1,y1)
    x1, y1 = 0.0, float(y[0])
    x2, y2 = float(n - 1), float(y[-1])
    denom = np.hypot(x2 - x1, y2 - y1) or 1.0
    dist = np.abs((x2 - x1) * (y1 - y) - (x1 - x) * (y2 - y1)) / denom
    cutoff = int(np.argmax(dist)) + 1
    ses = ranked[:cutoff]
    log.info("ROSE: %d super-enhancers called (inflection at rank %d of %d)",
             len(ses), cutoff, n)

    out = ROOT / "data/tracks/super_enhancer.bed"
    with open(out, "w") as fh:
        for i, (chrom, cs, ce, sig) in enumerate(ses, 1):
            fh.write(f"{chrom}\t{cs}\t{ce}\tSE_{i}\t{min(999, int(sig))}\t.\t{sig:.3f}\n")
    import hashlib
    entry = {
        "name": "SuperEnhancer_MM1S_ROSE",
        "display_name": "MM.1S super-enhancers (ROSE-style)",
        "path": "data/tracks/super_enhancer.bed",
        "url": "derived from " + primary.get("url", primary["name"]),
        "file_accession": "", "experiment": primary.get("experiment", ""),
        "assay": "derived: ROSE-style stitch of H3K27ac peaks", "target": "H3K27ac",
        "cell_type": "MM.1S", "assembly": "GRCh38", "category": "super_enhancer",
        "kind": "peaks", "output_type": f"derived ({method}; promoter TSS+/-2kb excluded; "
                                        "12.5kb stitching; inflection-point cutoff)",
        "md5": "",
    }
    import hashlib as _h
    entry["md5"] = _h.md5(out.read_bytes()).hexdigest()
    _upsert_track(reg, entry)
    (ROOT / "data/tracks/super_enhancer_method.json").write_text(json.dumps({
        "method": "ROSE-style", "stitch_gap_bp": 12500,
        "promoter_exclusion": "TSS +/- 2 kb (GENCODE)",
        "signal_method": method, "cutoff": "max distance from first-last chord",
        "cell_line": "MM.1S", "n_elements_ranked": len(ranked),
        "n_super_enhancers": len(ses),
        "source_experiment": primary.get("experiment", ""),
    }, indent=2))


# ---------------------------------------------------- GEO track lifting ---
def _detect_build(path: Path) -> str | None:
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i > 20:
                    break
                low = line.lower()
                if not line.startswith(("#", "track", "browser")):
                    if "\t" in line and line.split("\t")[0].startswith("chr"):
                        pass
                for pat, build in (("hg38", "hg38"), ("grch38", "hg38"),
                                   ("hg19", "hg19"), ("grch37", "hg19")):
                    if pat in low:
                        return build
    except Exception:  # noqa: BLE001
        return None
    return None


def register_geo_tracks(reg: dict) -> None:
    """Register GSE160335 control H3K27ac peaks; liftOver hg19 -> GRCh38."""
    geo_dir = ROOT / "data/tracks/geo"
    if not geo_dir.exists():
        log.info("no GEO track dir; skipping")
        return
    files = sorted(p for p in geo_dir.iterdir()
                   if p.suffix in (".bed", ".gz", ".txt", ".narrowPeak", ".broadPeak")
                   and "series_matrix" not in p.name and not p.name.startswith("GSE160335_RAW"))
    # plain .txt/.narrowPeak downloads from GEO suppl
    files += sorted(p for p in geo_dir.iterdir() if p.name.endswith((".narrowPeak", ".bed.txt")))
    files = [f for f in dict.fromkeys(files)]
    if not files:
        log.info("no GSE160335 control peak files found; skipping")
        return
    from pyliftover import LiftOver
    chain_path = ROOT / "data/reference/hg19ToHg38.over.chain.gz"
    lo = None
    lifted_dir = geo_dir / "lifted"
    lifted_dir.mkdir(exist_ok=True)
    n_ok = n_lifted = 0
    for f in files:
        build = _detect_build(f)
        out = lifted_dir / (f.name.replace(".gz", "") + ".hg38.bed")
        if build == "hg38":
            # normalize to plain bed in lifted dir for consistency
            out = lifted_dir / (f.name.replace(".gz", "") + ".hg38.bed")
            if not out.exists():
                opener = gzip.open if f.suffix == ".gz" else open
                with opener(f, "rt", errors="replace") as src, open(out, "w") as dst:
                    for line in src:
                        if line.startswith(("#", "track", "browser")) or not line.strip():
                            continue
                        dst.write(line)
            src_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160335/suppl/" + f.name
        else:
            if lo is None:
                if not chain_path.exists():
                    raise FileNotFoundError(f"liftOver chain missing: {chain_path}")
                lo = LiftOver(str(chain_path))
            if not out.exists():
                opener = gzip.open if f.suffix == ".gz" else open
                with opener(f, "rt", errors="replace") as src, open(out, "w") as dst:
                    for line in src:
                        if line.startswith(("#", "track", "browser")) or not line.strip():
                            continue
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 3:
                            continue
                        try:
                            s, e = int(parts[1]), int(parts[2])
                        except ValueError:
                            continue
                        c1 = lo.convert(parts[0], s)
                        c2 = lo.convert(parts[0], e - 1)
                        if not c1 or not c2 or c1[0] != c2[0]:
                            continue
                        nc = c1[0] if str(c1[0]).startswith("chr") else f"chr{c1[0]}"
                        ns, ne = sorted((int(c1[1]), int(c2[1]) + 1))
                        parts[0], parts[1], parts[2] = nc, str(ns), str(ne)
                        dst.write("\t".join(parts) + "\n")
                        n_lifted += 1
            src_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160335/suppl/" + f.name
        if out.exists() and out.stat().st_size > 100:
            import hashlib
            name = f"ControlH3K27ac_GSE160335_{f.stem.split('_')[-1][:12]}_{n_ok}"
            entry = {
                "name": name, "display_name": f"MM.1S control H3K27ac (GSE160335, {f.name[:24]})",
                "path": str(out.relative_to(ROOT)), "url": src_url,
                "file_accession": f"GSE160335/{f.name}", "experiment": "GSE160335",
                "assay": "H3K27ac ChIP-seq (untreated/control only)", "target": "H3K27ac",
                "cell_type": "MM.1S",
                "assembly": "GRCh38" if build == "hg38" else
                           "GRCh38 (liftOver from hg19 via hg19ToHg38.over.chain.gz)",
                "category": "h3k27ac", "kind": "peaks",
                "output_type": "processed peaks (control sample)",
                "md5": hashlib.md5(out.read_bytes()).hexdigest(),
            }
            _upsert_track(reg, entry)
            n_ok += 1
    log.info("GEO tracks registered: %d files (%d lifted intervals: %d)", n_ok, n_lifted)


# -------------------------------------------------------- track indexes ---
def build_track_indexes(cfg, reg: dict) -> None:
    from myelovar.io.bed import IntervalIndex, load_bed
    cache = ROOT / cfg.paths.cache_dir / "tracks"
    cache.mkdir(parents=True, exist_ok=True)
    for t in reg.get("tracks", []):
        if t.get("kind") != "peaks":
            continue
        dest = cache / f"{t['name']}.parquet"
        if dest.exists():
            continue
        df, _ = load_bed(ROOT / t["path"])
        df.to_parquet(dest, index=False)
        IntervalIndex(df, merge=True)
        log.info("track index built: %s (%d intervals)", t["name"], len(df))


def build_af_prewarm(cfg) -> None:
    """Precompute dbSNP-common AF for the default Mode A input VCF."""
    from myelovar.io.vcf import load_vcf
    from myelovar.reference import ReferenceBundle
    vcf = ROOT / "data/raw/giab/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    if not vcf.exists():
        log.warning("GIAB VCF missing; skipping AF prewarm")
        return
    t0 = time.time()
    df = load_vcf(vcf, autosomes_only=cfg.filter.autosomes_only)
    bundle = ReferenceBundle(cfg, root=ROOT)
    out = bundle.af(df)
    log.info("AF prewarm done in %.0fs: %d/%d matched",
             time.time() - t0, int(out["in_common_track"].sum()), len(out))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-af", action="store_true", help="skip slow AF prewarm")
    args = ap.parse_args()
    from myelovar.config import load_config
    cfg = load_config(ROOT / "config.yaml")

    registry_path = ROOT / cfg.paths.track_registry
    reg = _load_registry(registry_path)

    build_gencode(cfg)
    register_geo_tracks(reg)
    _save_registry(registry_path, reg)
    build_super_enhancer(cfg, reg)
    _save_registry(registry_path, reg)
    build_track_indexes(cfg, reg)
    build_abc(cfg)
    build_depmap(cfg)
    build_gwas_loci(cfg)
    if not args.skip_af:
        build_af_prewarm(cfg)
    log.info("all caches built")
    return 0


if __name__ == "__main__":
    sys.exit(main())
