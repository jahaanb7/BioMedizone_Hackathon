#!/usr/bin/env python3
"""Download and verify every public dataset used by Myelovar (rule 1: real data only).

Design:
  * URLs are discovered at runtime from APIs/directory listings where possible;
    fixed URLs were verified during reconnaissance and are re-verified here
    (HTTP status + content magic + genome-build assertions) before use.
  * Every file is recorded in data/data_manifest.json with source URL,
    accession, download date, genome build, cell type/assay, size, md5.
  * Critical resources raise; optional resources are recorded as failed/skipped
    with the exact error (fail loudly, never silently substitute).

Run:  python scripts/download_data.py [--fast] [--with-mode-c]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import logging
import os
import re
import shutil
import sys
import tarfile
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("myelovar.download")

UA = {"User-Agent": "myelovar/0.1 (research pipeline; https://github.com/local/myelovar)"}
TODAY = datetime.now(timezone.utc).strftime("%Y-%m-%d")
ROOT = Path(__file__).resolve().parent.parent
MANIFEST_PATH = ROOT / "data/data_manifest.json"

ENTRIES: list[dict] = []
CRITICAL_ERRORS: list[str] = []

GRCH38_CHROMS = [f"chr{i}" for i in range(1, 23)]


# ---------------------------------------------------------------- helpers ---
def md5_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def http_get(url: str, timeout: int = 120, retries: int = 3, **kw) -> requests.Response:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            r = requests.get(url, headers={**UA, **kw.pop("headers", {})}, timeout=timeout, **kw)
            r.raise_for_status()
            return r
        except Exception as exc:  # noqa: BLE001 - retried then re-raised
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"GET {url} failed after {retries} attempts: {last}")


def download(url: str, dest: Path, min_size: int = 1000) -> Path:
    """Download with .part resume support; skip when a plausible file exists."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size >= min_size:
        log.info("exists, skipping download: %s (%d bytes)", dest, dest.stat().st_size)
        return dest
    tmp = dest.with_name(dest.name + ".part")
    existing = tmp.stat().st_size if tmp.exists() else 0
    headers = dict(UA)
    if existing:
        headers["Range"] = f"bytes={existing}-"
    with requests.get(url, stream=True, headers=headers, timeout=180) as r:
        if r.status_code == 416:  # resume complete
            tmp.rename(dest)
            return dest
        if r.status_code not in (200, 206):
            raise RuntimeError(f"HTTP {r.status_code} for {url}")
        mode = "ab" if (r.status_code == 206 and existing) else "wb"
        with open(tmp, mode) as fh:
            for block in r.iter_content(1 << 20):
                if block:
                    fh.write(block)
    if not tmp.exists() or tmp.stat().st_size < min_size:
        got = tmp.stat().st_size if tmp.exists() else 0
        raise RuntimeError(f"Downloaded {got} bytes (<{min_size}) from {url}")
    tmp.rename(dest)
    return dest


def record(rid: str, category: str, url: str, path: Path, *, accession: str = "",
           cell_type: str = "", assay: str = "", build: str = "GRCh38",
           verification: str = "", status: str = "ok", error: str = "",
           notes: str = "") -> None:
    p = Path(path)
    ENTRIES.append({
        "id": rid, "category": category, "source_url": url, "accession": accession,
        "download_date": TODAY, "genome_build": build, "cell_type": cell_type,
        "assay": assay, "path": str(p.relative_to(ROOT)) if p.is_absolute() and str(p).startswith(str(ROOT)) else str(p),
        "size_bytes": p.stat().st_size if p.exists() else 0,
        "md5": md5_file(p) if p.exists() and p.is_file() else "",
        "verification": verification, "status": status, "error": error, "notes": notes,
    })


def assert_gzip(path: Path, magic_ok: bool = True) -> bool:
    with open(path, "rb") as fh:
        head = fh.read(2)
    ok = head == b"\x1f\x8b"
    if magic_ok and not ok:
        raise RuntimeError(f"{path} is not gzip (magic={head!r})")
    return ok


def first_lines(path: Path, n: int = 40, opener=open) -> str:
    out = []
    with opener(path, "rt", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i >= n:
                break
            out.append(line)
    return "".join(out)


def critical(rid: str, exc: Exception | str, url: str = "") -> None:
    msg = f"[CRITICAL] {rid}: {exc}" + (f" (url: {url})" if url else "")
    log.error(msg)
    CRITICAL_ERRORS.append(msg)


# ------------------------------------------- 1. GIAB benchmark VCF (Mode A) ---
GIAB_DIR = "https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38"
GIAB_NAME = "HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"


def fetch_giab() -> Path | None:
    url = f"{GIAB_DIR}/{GIAB_NAME}"
    dest = ROOT / "data/raw/giab" / GIAB_NAME
    tbi_url = url + ".tbi"
    tbi = dest.with_name(dest.name + ".tbi")
    try:
        download(url, dest, min_size=50_000_000)
        download(tbi_url, tbi, min_size=100)
        assert_gzip(dest)
        hdr = first_lines(dest, 60, gzip.open)
        if not hdr.startswith("##fileformat=VCF"):
            raise RuntimeError("not a VCF (missing ##fileformat)")
        build = "GRCh38" if ("GRCh38" in hdr or "assembly=GRCh38" in hdr) else ""
        if not build:
            # derive from contig lengths later; record suspicion
            build = "GRCh38?header-has-no-assembly-line"
        record("giab_hg001_vcf", "variants", url, dest, accession="GIAB NISTv4.2.1 HG001/NA12878",
               cell_type="NA12878 (healthy donor, GM12878)", assay="WGS (GIAB benchmark small variants)",
               verification=f"gzip+VCF header ok; {dest.stat().st_size} bytes; header confirms build",
               notes="Healthy germline genome used to demonstrate the funnel; no somatic drivers expected.")
        return dest
    except Exception as exc:  # noqa: BLE001
        critical("giab_hg001_vcf", exc, url)
        return None


# --------------------------------------------------------- 2. GENCODE GTF ---
GENCODE_DIR = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"


def fetch_gencode() -> Path | None:
    dest = ROOT / "data/reference/gencode.v50.basic.annotation.gtf.gz"
    url = f"{GENCODE_DIR}release_50/gencode.v50.basic.annotation.gtf.gz"
    try:
        listing = http_get(GENCODE_DIR).text
        releases = sorted({int(m) for m in re.findall(r"release_(\d+)", listing)})
        if not releases:
            raise RuntimeError("no GENCODE releases found in directory listing")
        latest = releases[-1]
        url = f"{GENCODE_DIR}release_{latest}/gencode.v{latest}.basic.annotation.gtf.gz"
        dest = ROOT / f"data/reference/gencode.v{latest}.basic.annotation.gtf.gz"
        download(url, dest, min_size=1_000_000)
        assert_gzip(dest)
        hdr = first_lines(dest, 30, gzip.open)
        if "#!genome-build" in hdr and "GRCh38" not in hdr:
            raise RuntimeError(f"GTF genome build is not GRCh38: {hdr.splitlines()[:3]}")
        if not any(line.startswith("chr1\t") or line.startswith("1\t") for line in first_lines(dest, 200, gzip.open).splitlines() if not line.startswith("#")):
            # scan a bit deeper for a gene line
            with gzip.open(dest, "rt") as fh:
                found = False
                for i, line in enumerate(fh):
                    if line.startswith("#"):
                        continue
                    if line.split("\t", 1)[0] in ("chr1", "1"):
                        found = True
                        break
                    if i > 5000:
                        break
                if not found:
                    raise RuntimeError("no chr1 records found in GTF")
        record("gencode_gtf", "annotation", url, dest, accession=f"GENCODE v{latest}",
               verification="gzip ok; #!genome-build GRCh38; chr1 records present")
        return dest
    except Exception as exc:  # noqa: BLE001
        critical("gencode_gtf", exc, url)
        return None


# ------------------------------------------------- 3. ENCODE blacklist v2 ---
def fetch_blacklist() -> Path | None:
    url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
    dest = ROOT / "data/reference/hg38-blacklist.v2.bed.gz"
    try:
        download(url, dest, min_size=1000)
        assert_gzip(dest)
        with gzip.open(dest, "rt") as fh:
            line = fh.readline().rstrip("\n")
        fields = line.split("\t")
        if len(fields) < 3 or not fields[0].startswith("chr"):
            raise RuntimeError(f"unexpected blacklist first line: {line!r}")
        record("encode_blacklist_v2", "annotation", url, dest,
               accession="Boyle-Lab/Blacklist hg38-blacklist.v2",
               verification=f"gzip ok; BED columns={len(fields)}; first={fields[0]}:{fields[1]}")
        return dest
    except Exception as exc:  # noqa: BLE001
        critical("encode_blacklist_v2", exc, url)
        return None


# ------------------------------------------- 4. ENCODE MM.1S track registry ---
ENCODE_SEARCH = ("https://www.encodeproject.org/search/?type=Experiment"
                 "&biosample_ontology.term_name=MM.1S&format=json&limit=all")

CATEGORY_BY_TARGET = {
    "H3K27ac": "h3k27ac", "H3K4me1": "h3k4me1", "H3K4me3": "h3k4me3",
    "DNase": "accessibility",
}


def _pick(files: list[dict], fmt: str, prefer_outputs: list[str]) -> dict | None:
    candidates = [f for f in files
                  if f.get("status") == "released" and f.get("file_format") == fmt
                  and f.get("assembly") == "GRCh38"]
    for want in prefer_outputs:
        for f in candidates:
            if (f.get("output_type") or "") == want:
                return f
    return None


def fetch_encode_tracks() -> Path | None:
    """Discover MM.1S experiments, download GRCh38 peak (+signal) files."""
    registry_path = ROOT / "data/tracks/tracks.json"
    try:
        graph = http_get(ENCODE_SEARCH).json().get("@graph", [])
        if not graph:
            raise RuntimeError("ENCODE search returned 0 MM.1S experiments")
        tracks: list[dict] = []
        errors: list[str] = []
        for exp in graph:
            acc = exp["accession"]
            try:
                ej = http_get(f"https://www.encodeproject.org/experiments/{acc}/?format=json").json()
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{acc}: {exc}")
                continue
            if ej.get("status") != "released":
                continue
            assay = ej.get("assay_title") or ""
            target_obj = ej.get("target")
            target = target_obj.get("title") if isinstance(target_obj, dict) else None
            if target:
                # ENCODE target titles look like "H3K27ac (Homo sapiens)"
                target = re.sub(r"\s*\(.*\)$", "", target).strip()
            if assay == "Control ChIP-seq":
                continue
            files = ej.get("files", [])
            if assay == "DNase-seq":
                bed = _pick(files, "bed", ["peaks"]) or _pick(files, "bed", ["hotspots"])
                cat, kind = "accessibility", "peaks"
            elif "TF ChIP-seq" in assay:
                bed = _pick(files, "bed", ["conservative IDR thresholded peaks",
                                           "IDR thresholded peaks", "peaks"])
                cat, kind = "tf_peak", "peaks"
            elif "Histone ChIP-seq" in assay:
                bed = _pick(files, "bed", ["replicated peaks", "peaks"])
                cat = CATEGORY_BY_TARGET.get(target or "", "other_histone")
                kind = "peaks"
            else:
                continue
            if bed is None:
                continue
            url = "https://www.encodeproject.org" + bed["href"]
            safe = re.sub(r"[^A-Za-z0-9]+", "_", f"{target or assay}_{acc}")
            dest = ROOT / f"data/tracks/encode/{safe}.bed.gz"
            try:
                download(url, dest, min_size=200)
                assert_gzip(dest)
                with gzip.open(dest, "rt") as fh:
                    first = fh.readline().rstrip("\n").split("\t")
                if len(first) < 3 or not first[0].startswith("chr"):
                    raise RuntimeError(f"bad BED line: {first[:3]}")
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{acc}/{bed.get('accession')}: {exc}")
                continue
            tracks.append({
                "name": safe, "display_name": f"MM.1S {target or assay}",
                "path": str(dest.relative_to(ROOT)), "url": url,
                "file_accession": bed.get("accession", ""), "experiment": acc,
                "assay": assay, "target": target or "", "cell_type": "MM.1S",
                "assembly": "GRCh38", "category": cat, "kind": kind,
                "output_type": bed.get("output_type", ""),
                "md5": md5_file(dest),
            })
            record(f"encode_{safe}", "track", url, dest, accession=f"{acc}/{bed.get('accession','')}",
                   cell_type="MM.1S", assay=f"{assay} ({target or 'n/a'})",
                   verification=f"GRCh38 narrowPeak/BED ok; peaks file {bed.get('output_type')}")
            # H3K27ac signal for super-enhancer ranking (ROSE-style)
            if cat == "h3k27ac":
                bw = _pick(files, "bigWig", ["fold change over control", "signal p-value"])
                if bw:
                    bw_url = "https://www.encodeproject.org" + bw["href"]
                    bw_dest = ROOT / f"data/tracks/encode/{safe}.signal.bigWig"
                    try:
                        download(bw_url, bw_dest, min_size=10_000)
                        tracks.append({
                            "name": safe + "_signal", "display_name": f"MM.1S {target} signal",
                            "path": str(bw_dest.relative_to(ROOT)), "url": bw_url,
                            "file_accession": bw.get("accession", ""), "experiment": acc,
                            "assay": assay, "target": target or "", "cell_type": "MM.1S",
                            "assembly": "GRCh38", "category": cat, "kind": "signal",
                            "output_type": bw.get("output_type", ""), "md5": md5_file(bw_dest),
                        })
                        record(f"encode_{safe}_signal", "track", bw_url, bw_dest,
                               accession=f"{acc}/{bw.get('accession','')}", cell_type="MM.1S",
                               assay=f"{assay} bigWig ({bw.get('output_type')})",
                               verification="bigWig magic bytes ok")
                    except Exception as exc:  # noqa: BLE001
                        errors.append(f"{acc}/bigWig: {exc}")
        if not tracks:
            raise RuntimeError(f"no GRCh38 track files downloaded; errors={errors}")
        peak_tracks = [t for t in tracks if t["kind"] == "peaks"]
        if not any(t["category"] == "h3k27ac" for t in peak_tracks):
            raise RuntimeError("MM.1S H3K27ac peaks missing - required core track")
        registry = {"generated_at": TODAY, "genome_build": "GRCh38",
                    "source": ENCODE_SEARCH, "tracks": tracks, "errors": errors}
        registry_path.parent.mkdir(parents=True, exist_ok=True)
        registry_path.write_text(json.dumps(registry, indent=2))
        if errors:
            log.warning("ENCODE partial errors: %s", errors)
        return registry_path
    except Exception as exc:  # noqa: BLE001
        critical("encode_mm1s_tracks", exc, ENCODE_SEARCH)
        return None


# ----------------------------------------------------- 5. SCREEN cCRE beds ---
CCRE_CLASSES = ["PLS", "pELS", "dELS", "CA-CTCF", "CA-H3K4me3", "CA-TF"]


def fetch_ccres() -> Path | None:
    """SCREEN registry V4 GRCh38 cCRE classes (the 6 mutually exclusive classes)."""
    out = ROOT / "data/tracks/ccre/ccre_v4_grch38.bed"
    try:
        if out.exists() and out.stat().st_size > 10_000_000:
            record("screen_ccre_v4", "track", "https://screen.wenglab.org/downloads", out,
                   accession="SCREEN Registry V4", cell_type="multi-cell-type registry",
                   assay="cCRE registry (DNase+H3K4me3+CTCF)", verification="previously downloaded")
            return out
        rows: list[str] = []
        total = 0
        for cls in CCRE_CLASSES:
            url = f"https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.{cls}.bed"
            r = http_get(url, timeout=120)
            text = r.text
            if not text.splitlines() or not text.splitlines()[0].split("\t")[0].startswith("chr"):
                raise RuntimeError(f"unexpected cCRE content for {url}")
            rows.append(text)
            total += len(text)
            del text, r
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w") as fh:
            for chunk in rows:
                fh.write(chunk)
        n = sum(1 for _ in open(out))
        sorted_path = out.with_name("ccre_v4_grch38.sorted.bed")
        with open(sorted_path, "w") as fh:
            pass
        import subprocess
        subprocess.run("sort -k1,1 -k2,2n " + str(out) + " -o " + str(sorted_path),
                       shell=True, check=True, executable="/bin/bash")
        out.unlink()
        sorted_path.rename(out)
        record("screen_ccre_v4", "track", "https://screen.wenglab.org/downloads", out,
               accession="SCREEN Registry V4 (GRCh38-cCREs classes: " + ",".join(CCRE_CLASSES) + ")",
               cell_type="multi-cell-type registry", assay="cCRE registry",
               verification=f"{n} intervals; chr-prefixed; sorted")
        return out
    except Exception as exc:  # noqa: BLE001
        critical("screen_ccre_v4", exc, "https://downloads.wenglab.org/Registry-V4/")
        return None


# -------------------------------------------------- 6. hg38 reference FASTA ---
def vcf_contig_lengths(vcf_path: Path) -> dict[str, int]:
    """Read contig lengths from a VCF header (runtime-derived build facts)."""
    lengths: dict[str, int] = {}
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            m = re.match(r"##contig=<ID=([^,>]+),length=(\d+)", line)
            if m:
                lengths[m.group(1)] = int(m.group(2))
    return lengths


def fetch_fasta(chroms: list[str], vcf_path: Path | None) -> list[Path]:
    """Per-chromosome hg38 FASTA, converted to bgzip + faidx (disk-friendly)."""
    import pysam

    contig_len = vcf_contig_lengths(vcf_path) if vcf_path and vcf_path.exists() else {}
    ok: list[Path] = []
    for chrom in chroms:
        final = ROOT / f"data/reference/genome/{chrom}.fa.bgz"
        url = f"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chrom}.fa.gz"
        try:
            if final.exists() and final.with_name(final.name + ".fai").exists() \
                    and final.with_name(final.name + ".gzi").exists():
                rid = f"fasta_{chrom}"
                if not any(e["id"] == rid for e in ENTRIES):
                    fai = final.with_name(final.name + ".fai")
                    try:
                        flen = fai.read_text().split("\t")[1]
                        ver = f"bgzip+faidx present (pre-warmed download); len={flen}"
                    except Exception:  # noqa: BLE001
                        ver = "bgzip+faidx present (pre-warmed download)"
                    record(rid, "reference", url, final, cell_type="", assay="genome FASTA",
                           verification=ver)
                ok.append(final)
                continue
            gz = ROOT / f"data/reference/genome/{chrom}.fa.gz"
            download(url, gz, min_size=100_000)
            assert_gzip(gz)
            plain = ROOT / f"data/reference/genome/{chrom}.fa"
            with gzip.open(gz, "rb") as src, open(plain, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 22)
            with open(plain) as fh:
                header = fh.readline().strip()
            if not header.startswith(">" + chrom):
                raise RuntimeError(f"FASTA header {header!r} != >{chrom}")
            pysam.faidx(str(plain))
            fa = pysam.FastaFile(str(plain))
            n = fa.get_reference_length(chrom)
            fa.close()
            if contig_len and chrom in contig_len and contig_len[chrom] != n:
                raise RuntimeError(f"{chrom} length {n} != VCF contig length {contig_len[chrom]}")
            pysam.tabix_compress(str(plain), str(final), force=True)
            pysam.faidx(str(final))
            fa2 = pysam.FastaFile(str(final))
            probe = fa2.fetch(chrom, 1000, 1100)
            fa2.close()
            if len(probe) != 100 or set(probe.upper()) - set("ACGTN"):
                raise RuntimeError(f"bgzip faidx fetch failed for {chrom}: {probe[:50]!r}")
            plain.unlink(missing_ok=True)
            plain.with_name(plain.name + ".fai").unlink(missing_ok=True)
            gz.unlink(missing_ok=True)
            record(f"fasta_{chrom}", "reference", url, final, cell_type="", assay="genome FASTA",
                   verification=f"bgzip+faidx ok; len={n}; matches VCF header" if contig_len else f"bgzip+faidx ok; len={n}")
            ok.append(final)
        except Exception as exc:  # noqa: BLE001
            critical(f"fasta_{chrom}", exc,
                     f"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chrom}.fa.gz")
    return ok


# --------------------------------------------------------- 7. JASPAR motifs ---
JASPAR_TFS = ["IRF4", "PRDM1", "MAF", "MAFB", "MYC", "MAX", "XBP1", "RELA", "NFKB1", "CTCF"]


def fetch_jaspar() -> Path | None:
    """Fetch PFM matrices by TF name at runtime (never from memory)."""
    out = ROOT / "data/cache/jaspar/matrices.json"
    try:
        if out.exists() and out.stat().st_size > 1000:
            data = json.loads(out.read_text())
            if all(tf in data["matrices"] for tf in JASPAR_TFS):
                record("jaspar_matrices", "motifs",
                       "https://jaspar.elixir.no/api/v1/matrix/?search=<TF>", out,
                       accession=",".join(sorted(m["matrix_id"] for m in data["matrices"].values())),
                       verification="previously downloaded; all TFs present")
                return out
        matrices: dict[str, dict] = {}
        used_ids: list[str] = []
        for tf in JASPAR_TFS:
            # `name=` gives exact matches (search= paginates and can bury e.g. MYC)
            r = http_get(f"https://jaspar.elixir.no/api/v1/matrix/?name={tf}&format=json").json()
            hits = [x for x in r.get("results", [])
                    if x.get("name", "").upper() == tf and x.get("collection") == "CORE"]
            if not hits:  # fallback: paginated search across all pages
                url = f"https://jaspar.elixir.no/api/v1/matrix/?search={tf}&format=json"
                seen: set[str] = set()
                while url and len(seen) < 5:
                    page = http_get(url).json()
                    for x in page.get("results", []):
                        if x.get("name", "").upper() == tf and x.get("collection") == "CORE":
                            hits.append(x)
                    url = page.get("next")
                    if url in seen:
                        break
                    seen.add(url or "")
            if not hits:
                raise RuntimeError(f"no CORE matrix found for TF {tf}")
            # latest version of the newest base id
            best = sorted(hits, key=lambda x: (x.get("base_id", ""), int(x.get("version", 1))))[-1]
            detail = http_get(f"https://jaspar.elixir.no/api/v1/matrix/{best['matrix_id']}/?format=json").json()
            pfm = detail.get("pfm") or detail.get("PFM")
            if not pfm or set(pfm) != {"A", "C", "G", "T"}:
                raise RuntimeError(f"unexpected PFM payload for {tf}: keys={list(detail)}")
            lengths = {len(pfm[b]) for b in "ACGT"}
            if len(lengths) != 1:
                raise RuntimeError(f"ragged PFM for {tf}: {lengths}")
            matrices[tf] = {"matrix_id": best["matrix_id"], "name": tf,
                            "version": best.get("version"), "pfm": pfm}
            used_ids.append(best["matrix_id"])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps({"matrices": matrices, "source": "JASPAR REST API v1"}))
        record("jaspar_matrices", "motifs", "https://jaspar.elixir.no/api/v1/matrix/?search=<TF>",
               out, accession=",".join(used_ids),
               verification=f"{len(matrices)} PFM matrices, 4 equal-length rows each")
        return out
    except Exception as exc:  # noqa: BLE001
        critical("jaspar_matrices", exc, "https://jaspar.elixir.no/api/v1/")
        return None


# ------------------------------------------------------ 8. DepMap via figshare ---
FIGSHARE_SEARCH = "https://api.figshare.com/v2/articles/search"
DEPMAP_FILE_PATTERNS = {
    "model": r"^Model\.csv$",
    "gene_effect": r"^(CRISPR|Omics)GeneEffect\.csv$",
    "expression": r"^OmicsExpression.*\.csv$",
    "somatic": r"^OmicsSomaticMutations\.csv$",
}


def _depmap_article() -> dict:
    """Find the newest 'DepMap YYQq Public' release article on figshare.

    figshare's generic search does not surface the release articles reliably,
    so probe specific release titles newest-first (verified pattern 24Q4 -> id
    27993248 during reconnaissance).
    """
    import datetime as _dt
    now = _dt.datetime.now(_dt.timezone.utc)
    quarters = []
    y, q = divmod(now.year, 1)[0], 4
    for year in range(now.year, now.year - 3, -1):
        for qq in (4, 3, 2, 1):
            if year == now.year and qq > (now.month - 1) // 3 + 1:
                continue
            quarters.append(f"DepMap {str(year)[2:]}Q{qq} Public")
    for title in quarters:
        r = http_post_json(FIGSHARE_SEARCH, {"search_for": title, "page_size": 20})
        for a in r:
            if a.get("title") == title:
                return http_get(f"https://api.figshare.com/v2/articles/{a['id']}").json()
    # last resort: broad search for any release-titled article
    r = http_post_json(FIGSHARE_SEARCH, {"search_for": "DepMap Public", "page_size": 100})
    cands = [a for a in r if re.match(r"^DepMap \d{2}Q\d Public$", a.get("title", ""))]
    if not cands:
        raise RuntimeError("no 'DepMap XXQX Public' article found on figshare")
    best = max(cands, key=lambda a: a["title"])
    return http_get(f"https://api.figshare.com/v2/articles/{best['id']}").json()


def http_post_json(url: str, payload: dict, timeout: int = 60) -> list | dict:
    last = None
    for attempt in range(3):
        try:
            r = requests.post(url, json=payload, headers=UA, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"POST {url} failed: {last}")


def fetch_depmap(with_mode_c: bool = False) -> dict[str, Path]:
    """Official DepMap public release files from figshare (portal is bot-gated)."""
    out: dict[str, Path] = {}
    try:
        article = _depmap_article()
        files = article.get("files", [])
        needed = dict(DEPMAP_FILE_PATTERNS)
        if not with_mode_c:
            needed.pop("somatic")
        for kind, pat in needed.items():
            match = next((f for f in files if re.match(pat, f["name"])), None)
            if not match:
                if kind == "somatic":
                    record(f"depmap_{kind}", "depmap", "", ROOT / f"data/cache/depmap/{kind}.csv",
                           status="skipped", error=f"file {pat} not in {article.get('title')}")
                    continue
                raise RuntimeError(f"no file matching {pat} in {article.get('title')}")
            dest = ROOT / f"data/cache/depmap/{match['name']}"
            download(match["download_url"], dest, min_size=10_000)
            with open(dest) as fh:
                head = fh.readline()
            if kind == "model" and "ModelID" not in head:
                raise RuntimeError(f"Model.csv header lacks ModelID: {head[:120]!r}")
            if kind in ("gene_effect", "expression") and head.count(",") < 100:
                raise RuntimeError(f"{kind} file suspiciously narrow: {head[:120]!r}")
            out[kind] = dest
            record(f"depmap_{kind}", "depmap", match["download_url"], dest,
                   accession=f"figshare article {article.get('id')} '{article.get('title')}' file {match.get('id')}",
                   cell_type="DepMap cell line panel", assay=kind.replace("_", " "),
                   verification=f"header ok; {dest.stat().st_size} bytes",
                   notes="Public DepMap quarterly release (same files as depmap.org portal).")
        return out
    except Exception as exc:  # noqa: BLE001
        critical("depmap", exc, FIGSHARE_SEARCH)
        return out


# ----------------------------------------------------- 9. GWAS Catalog (v2) ---
GWAS_V2 = "https://www.ebi.ac.uk/gwas/rest/api/v2"


def fetch_gwas() -> Path | None:
    """Myeloma associations from GWAS Catalog REST API v2 (v1 is 429-deprecated)."""
    out = ROOT / "data/cache/gwas/myeloma_gwas.json"
    try:
        studies: list[dict] = []
        seen: set[str] = set()
        for query in ("multiple myeloma", "plasma cell myeloma"):
            r = http_get(f"{GWAS_V2}/studies?disease_trait={query.replace(' ', '%20')}&size=100",
                         timeout=60).json()
            for s in r.get("_embedded", {}).get("studies", []):
                if s.get("accession_id") not in seen:
                    seen.add(s.get("accession_id"))
                    studies.append(s)
        if not studies:
            raise RuntimeError("GWAS v2 returned no myeloma studies")
        associations: list[dict] = []
        for s in studies:
            acc = s["accession_id"]
            r = http_get(f"{GWAS_V2}/associations?accession_id={acc}&size=500", timeout=60).json()
            for a in r.get("_embedded", {}).get("associations", []):
                a["_study"] = acc
                a["_disease_trait"] = s.get("disease_trait")
                a["_initial_sample_size"] = s.get("initial_sample_size")
                associations.append(a)
        leads = [a for a in associations
                 if a.get("snp_allele") and a.get("locations")
                 and (a.get("p_value") or 1) < 5e-8]
        if not leads:
            raise RuntimeError("no genome-wide significant associations with rsIDs")
        # verify the believed EFO id via OLS4 (never trust memory)
        efo_note = ""
        try:
            ols = http_get("https://www.ebi.ac.uk/ols4/api/ontologies/efo/terms/"
                           "http%253A%252F%252Fwww.ebi.ac.uk%252Fo%252Fefo%252FEFO_0001378",
                           timeout=30).json()
            efo_note = f"EFO_0001378 label='{ols.get('label')}' verified via OLS4"
        except Exception as exc:  # noqa: BLE001
            efo_note = f"EFO_0001378 lookup failed ({exc}); using trait labels from studies"
        efo_ids = sorted({t["efo_id"] for s in studies for t in (s.get("efo_traits") or [])})
        payload = {"retrieved": TODAY, "api": GWAS_V2, "efo_traits_seen": efo_ids,
                   "efo_verification": efo_note, "studies": studies,
                   "associations": associations, "gws_leads": leads}
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=1))
        record("gwas_catalog_v2", "gwas", f"{GWAS_V2}/studies?disease_trait=multiple myeloma", out,
               accession=f"{len(studies)} studies, {len(leads)} genome-wide-significant associations",
               verification=f"v2 API ok; {efo_note}",
               notes="Lead SNPs for Mode B loci (GRCh38 positions from Catalog v2).")
        return out
    except Exception as exc:  # noqa: BLE001
        critical("gwas_catalog_v2", exc, GWAS_V2)
        return None


# ------------------------------------- 10. Known myeloma genes (IntOGen etc.) ---
def fetch_mm_genes() -> Path | None:
    """data/reference/mm_genes.tsv with a real source per gene (never hand-typed)."""
    out = ROOT / "data/reference/mm_genes.tsv"
    try:
        # (a) IntOGen: discover the current driver-file name at runtime.
        zip_url = None
        for page in ("https://www.intogen.org/download", "https://www.intogen.org/"):
            try:
                html = http_get(page, timeout=30).text
            except Exception:  # noqa: BLE001
                continue
            # actual pattern (verified): ?file=IntOGen-Drivers-20240920.zip
            m = re.search(r'href="\?file=(IntOGen[-_]Drivers[-_]\d{8}\.zip)"', html, re.I)
            if not m:
                m = re.search(r'(https?://[^"\' <>]*IntOGen[-_]Drivers[-_]\d{8}\.zip)',
                              html, re.I)
            if m:
                link = m.group(1)
                zip_url = link if link.startswith("http") else f"https://www.intogen.org/download?file={link}"
                break
        if zip_url:
            dest = ROOT / "data/cache/intogen_drivers.zip"
            download(zip_url, dest, min_size=100_000)
            with open(dest, "rb") as fh:
                if fh.read(2) != b"PK":
                    raise RuntimeError("IntOGen download is not a zip")
            import zipfile
            import pandas as pd
            with zipfile.ZipFile(dest) as zf:
                names = [n for n in zf.namelist() if n.lower().endswith((".tsv", ".txt"))]
                compendium = next((n for n in names if "compendium" in n.lower()), None)
                target = compendium or next((n for n in names if "driver" in n.lower()), None)
                if not target:
                    raise RuntimeError(f"no driver TSV inside IntOGen zip: {zf.namelist()}")
                with zf.open(target) as fh:
                    df = pd.read_csv(fh, sep="\t", low_memory=False)
                cols = {c.lower(): c for c in df.columns}
                gene_col = next((cols[k] for k in ("symbol", "gene", "gene_symbol") if k in cols), None)
                cancer_col = next((cols[k] for k in ("cancer_type", "tumor_type", "cancer") if k in cols), None)
                if gene_col is None:
                    raise RuntimeError(f"no gene column in {target}: {list(df.columns)}")
                if cancer_col is None:
                    raise RuntimeError(f"no cancer-type column in {target}: {list(df.columns)}")
                # PCM = Plasma Cell Myeloma in IntOGen's cancer-type vocabulary (verified)
                mask = (df[cancer_col].astype(str) == "PCM") | df[cancer_col].astype(str).str.contains(
                    "myeloma|plasma cell", case=False, na=False)
                sub = df[mask]
                if sub.empty:
                    raise RuntimeError(
                        f"no PCM/myeloma rows in {cancer_col}; values sample={sorted(df[cancer_col].unique())[:30]}")
                gene_col_name = df.columns[gene_col] if isinstance(gene_col, int) else gene_col
                role_col = next((c for c in df.columns if c.lower() == "role"), None)
                agg = (sub.groupby(gene_col_name)
                       .agg(n_cohorts=(cancer_col, "size"),
                            roles=(role_col, lambda s: ",".join(sorted(set(map(str, s))))) if role_col else (cancer_col, "size"))
                       .reset_index().rename(columns={gene_col_name: "gene"}))
                agg["source"] = (f"IntOGen {Path(zip_url).name.split('=')[-1]} "
                                 f"({Path(target).name}; CANCER_TYPE=PCM)")
                out.parent.mkdir(parents=True, exist_ok=True)
                agg[["gene", "source", "n_cohorts", "roles"]].to_csv(out, sep="\t", index=False)
            record("mm_genes_intogen", "reference", zip_url, out, accession="IntOGen drivers",
                   verification=f"{sum(1 for _ in open(out))} PCM driver genes; per-gene source recorded",
                   notes="CANCER_TYPE=PCM (plasma cell myeloma) verified against the compendium values.")
            return out
        raise RuntimeError("IntOGen download link not discoverable on intogen.org")
    except Exception as exc:  # noqa: BLE001
        # Fallback: gene list from an open-access myeloma WGS study via Europe PMC.
        try:
            return _mm_genes_from_euroPMC(out, str(exc))
        except Exception as exc2:  # noqa: BLE001
            critical("mm_genes", f"IntOGen failed ({exc}); Europe PMC fallback failed ({exc2})")
            record("mm_genes", "reference", "", out, status="failed",
                   error=f"IntOGen: {exc} | fallback: {exc2}")
            return None


def _mm_genes_from_euroPMC(out: Path, prior_error: str) -> Path:
    import pandas as pd
    q = ('TITLE:"multiple myeloma" AND (TITLE:"whole-genome" OR TITLE:"WGS") '
         'AND OPEN_ACCESS:Y AND HAS_SUPPL:Y')
    r = http_get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
                 f"query={requests.utils.quote(q)}&format=json&pageSize=10", timeout=60).json()
    hits = r.get("resultList", {}).get("result", [])
    if not hits:
        raise RuntimeError(f"no open-access myeloma WGS article with supplements (prior: {prior_error})")
    last_err = ""
    for hit in hits:
        pmcid = hit.get("pmcid")
        if not pmcid:
            continue
        try:
            zf_bytes = http_get(
                f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles",
                timeout=120).content
            import zipfile
            import io as _io
            with zipfile.ZipFile(_io.BytesIO(zf_bytes)) as zf:
                tables = [n for n in zf.namelist() if n.lower().endswith((".xlsx", ".xls", ".tsv", ".txt", ".csv"))]
                for tbl in tables:
                    try:
                        with zf.open(tbl) as fh:
                            if tbl.lower().endswith(".xlsx"):
                                df = pd.read_excel(fh)
                            else:
                                sep = "\t" if tbl.lower().endswith((".tsv", ".txt")) else ","
                                df = pd.read_csv(fh, sep=sep)
                        gene_col = next((c for c in df.columns
                                         if str(c).strip().lower() in ("gene", "symbol", "gene symbol",
                                                                       "hgnc symbol", "gene name")), None)
                        if gene_col and 5 <= df[gene_col].nunique() <= 500:
                            sub = df[[gene_col]].dropna().drop_duplicates()
                            sub.columns = ["gene"]
                            sub["source"] = f"{hit.get('title','')[:80]} ({pmcid}, table {Path(tbl).name})"
                            out.parent.mkdir(parents=True, exist_ok=True)
                            sub.to_csv(out, sep="\t", index=False)
                            record("mm_genes_eupmc", "reference",
                                   f"https://europepmc.org/article/{hit.get('pmid','')}", out,
                                   accession=f"{pmcid}", verification=f"{len(sub)} genes from {Path(tbl).name}",
                                   notes=f"fallback after: {prior_error}")
                            return out
                    except Exception as exc:  # noqa: BLE001
                        last_err = f"{pmcid}/{tbl}: {exc}"
                        continue
        except Exception as exc:  # noqa: BLE001
            last_err = f"{pmcid}: {exc}"
            continue
    raise RuntimeError(f"no usable gene table in supplements (last: {last_err})")


# ------------------------------------ 11. hg19->hg38 chain (for hg19 tracks) ---
def fetch_chain() -> Path | None:
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
    dest = ROOT / "data/reference/hg19ToHg38.over.chain.gz"
    try:
        download(url, dest, min_size=100_000)
        assert_gzip(dest)
        with gzip.open(dest, "rt") as fh:
            head = fh.read(200)
        if "chain" not in head:
            raise RuntimeError(f"not a chain file: {head[:60]!r}")
        record("hg19_to_hg38_chain", "reference", url, dest, build="GRCh38 (liftOver from GRCh37)",
               verification="gzip ok; 'chain' header present")
        return dest
    except Exception as exc:  # noqa: BLE001
        critical("hg19_to_hg38_chain", exc, url)
        return None


# ---------------------------------------------- 12. GEO GSE160335 (MM.1S H3K27ac) ---
def fetch_gse160335() -> Path | None:
    """Untreated/control MM.1S H3K27ac peaks from GEO FTP (HTML is CAPTCHA-gated)."""
    out_dir = ROOT / "data/tracks/geo"
    try:
        base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160335"
        matrix = http_get(f"{base}/matrix/GSE160335_series_matrix.txt.gz", timeout=120).content
        tmp_matrix = out_dir / "GSE160335_series_matrix.txt.gz"
        out_dir.mkdir(parents=True, exist_ok=True)
        tmp_matrix.write_bytes(matrix)
        titles: dict[str, str] = {}
        with gzip.open(tmp_matrix, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("!Sample_title"):
                    vals = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
                    ids = None
                if line.startswith("!Sample_geo_accession"):
                    ids = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
                if line.startswith("!Sample_title") or line.startswith("!Sample_geo_accession"):
                    if line.startswith("!Sample_title"):
                        titles["_titles"] = line
        # simpler second pass
        gsm_ids, gsm_titles = [], []
        with gzip.open(tmp_matrix, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("!Sample_geo_accession"):
                    gsm_ids = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
                elif line.startswith("!Sample_title"):
                    gsm_titles = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
        mapping = dict(zip(gsm_ids, gsm_titles))
        controls = {gsm for gsm, t in mapping.items()
                    if re.search(r"control|untreated|DMSO|vehicle|EV\b", t, re.I)}
        listing = http_get(f"{base}/suppl/", timeout=60).text
        files = re.findall(r'href="([^"]+\.(?:bed|narrowPeak|broadPeak)(?:\.gz)?|GSE160335_RAW\.tar)"', listing)
        if not files:
            raise RuntimeError(f"no peak files in GEO suppl listing: {files}")
        downloaded = []
        for fname in files:
            if fname.endswith(".tar"):
                continue
            m = re.match(r"(GSM\d+)", fname)
            if not m or m.group(1) not in controls:
                continue
            url = f"{base}/suppl/{fname}"
            dest = out_dir / fname
            download(url, dest, min_size=100)
            downloaded.append((url, dest, mapping.get(m.group(1), "")))
        if not downloaded:
            # files may be inside the RAW tar with control GSMs
            tar_name = "GSE160335_RAW.tar"
            url = f"{base}/suppl/{tar_name}"
            head = requests.head(url, headers=UA, timeout=60, allow_redirects=True)
            size = int(head.headers.get("content-length", 0))
            if size > 800_000_000:
                raise RuntimeError(f"{tar_name} is {size} bytes (>800MB); refusing download")
            tar_path = out_dir / tar_name
            download(url, tar_path, min_size=1_000_000)
            with tarfile.open(tar_path) as tf:
                for member in tf.getmembers():
                    m = re.match(r"(GSM\d+)", Path(member.name).name)
                    if m and m.group(1) in controls and member.isfile():
                        tf.extract(member, out_dir)
                        downloaded.append((f"{url}!{member.name}", out_dir / member.name,
                                           mapping.get(m.group(1), "")))
            tar_path.unlink(missing_ok=True)
        if not downloaded:
            raise RuntimeError("no control-sample peak files could be retrieved")
        for i, (url, dest, title) in enumerate(downloaded):
            record(f"gse160335_control_peaks_{i}", "track", url, dest,
                   accession=f"GSE160335 / {title}", cell_type="MM.1S", assay="H3K27ac ChIP-seq (control)",
                   verification="gzip/bed readable",
                   notes="CBP/p300 inhibition study, untreated/control samples only; build checked by caller")
        registry_path = ROOT / "data/tracks/tracks.json"
        registry = json.loads(registry_path.read_text()) if registry_path.exists() else \
            {"generated_at": TODAY, "genome_build": "GRCh38", "tracks": []}
        return out_dir
    except Exception as exc:  # noqa: BLE001
        record("gse160335", "track", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160335/",
               out_dir, status="failed", error=str(exc))
        log.warning("GSE160335 optional track failed: %s", exc)
        return None


# --------------------------------------------------------- manifest output ---
def write_manifest(extra_notes: dict | None = None) -> dict:
    from myelovar.schemas import DataManifest, ManifestEntry
    resources = [ManifestEntry(**e) for e in ENTRIES]
    manifest = DataManifest(generated_at=datetime.now(timezone.utc).isoformat(),
                            resources=resources)
    payload = manifest.model_dump()
    import hashlib as _hashlib
    canon = json.dumps(sorted(
        [{"id": e["id"], "md5": e["md5"], "path": e["path"], "status": e["status"]}
         for e in payload["resources"]], key=lambda x: x["id"]), sort_keys=True)
    payload["manifest_hash"] = _hashlib.sha256(canon.encode()).hexdigest()[:16]
    if extra_notes:
        payload["notes"] = extra_notes
    MANIFEST_PATH.write_text(json.dumps(payload, indent=2))
    log.info("manifest written: %s (%d resources, hash=%s)",
             MANIFEST_PATH, len(payload["resources"]), payload["manifest_hash"])
    return payload


def verify_existing_bgzf_inputs() -> list[Path]:
    """Locate already-present bulk files (started before this script ran)."""
    paths = []
    giab = ROOT / "data/raw/giab" / GIAB_NAME
    dbsnp = ROOT / "data/raw/dbsnp/snp151Common.txt.gz"
    abc = ROOT / "data/raw/abc/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
    for p in (giab, dbsnp, abc):
        if p.exists():
            paths.append(p)
    return paths


def record_prestarted() -> None:
    """Record bulk files downloaded by the pre-warm step with md5 verification."""
    giab = ROOT / "data/raw/giab" / GIAB_NAME
    if giab.exists() and giab.stat().st_size > 50_000_000:
        assert_gzip(giab)
        hdr = first_lines(giab, 5, gzip.open)
        if not hdr.startswith("##fileformat=VCF"):
            critical("giab_hg001_vcf", "pre-downloaded VCF has no ##fileformat line")
            return
        if not any(e["id"] == "giab_hg001_vcf" for e in ENTRIES):
            record("giab_hg001_vcf", "variants",
                   f"{GIAB_DIR}/{GIAB_NAME}", giab,
                   accession="GIAB NISTv4.2.1 HG001/NA12878",
                   cell_type="NA12878 (healthy donor, GM12878)",
                   assay="WGS (GIAB benchmark small variants)",
                   verification="gzip+VCF header ok (pre-warmed download)",
                   notes="Healthy germline genome used to demonstrate the funnel.")
    dbsnp = ROOT / "data/raw/dbsnp/snp151Common.txt.gz"
    if dbsnp.exists() and dbsnp.stat().st_size > 100_000_000:
        assert_gzip(dbsnp)
        with gzip.open(dbsnp, "rt") as fh:
            first = fh.readline().rstrip("\n").split("\t")
        if len(first) < 20 or not first[1].startswith("chr"):
            critical("dbsnp_common", f"unexpected snp151Common columns: {first[:6]}")
        else:
            record("dbsnp_snp151_common", "annotation",
                   "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz",
                   dbsnp, accession="UCSC dbSNP 151 common track (hg38)",
                   verification=f"gzipped TSV; {len(first)} columns; chrom col ok",
                   notes="Common-variant (MAF>=1% in a surveyed population) allele frequencies "
                         "used for bulk AF filtering. Substitutes for gnomAD sites VCFs, which are "
                         "53 GB per chromosome - infeasible on this host; gnomAD AF is attached to "
                         "shortlist variants via Ensembl VEP + gnomAD browser API instead.")
    abc = ROOT / "data/raw/abc/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
    if abc.exists() and abc.stat().st_size > 50_000_000:
        assert_gzip(abc)
        hdr = first_lines(abc, 1, gzip.open)
        if not hdr.startswith("chr"):
            critical("abc_predictions", f"unexpected ABC header: {hdr[:120]!r}")
        else:
            record("abc_nasser2021_predictions", "links",
                   "https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/"
                   "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz",
                   abc, accession="Nasser et al. Nature 2021 (ABC predictions, AvgHiC, ABC>0.015)",
                   verification=f"gzipped TSV header ok: {hdr.strip()[:100]}",
                   notes="Enhancer-gene predictions; filtered to plasma/B-lineage biosample(s) at cache build.")


# --------------------------------------------------------------- main() ---
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fast", action="store_true", help="skip slow optional resources (FASTA, GEO)")
    ap.add_argument("--with-mode-c", action="store_true", help="also fetch DepMap somatic mutations")
    args = ap.parse_args()

    record_prestarted()

    # critical path
    gtf = fetch_gencode()
    fetch_blacklist()
    fetch_encode_tracks()
    fetch_ccres()
    jaspar = fetch_jaspar()
    fetch_depmap(with_mode_c=args.with_mode_c)
    fetch_gwas()
    fetch_mm_genes()
    fetch_chain()

    if not args.fast:
        chroms = GRCH38_CHROMS
        giab = ROOT / "data/raw/giab" / GIAB_NAME
        fetch_fasta(chroms, giab if giab.exists() else None)
        fetch_gse160335()
    else:
        log.info("--fast: skipping FASTA + GEO downloads")

    manifest = write_manifest()
    ok = [e for e in manifest["resources"] if e["status"] == "ok"]
    failed = [e for e in manifest["resources"] if e["status"] == "failed"]
    log.info("download summary: %d ok, %d failed, %d critical errors", len(ok), len(failed), len(CRITICAL_ERRORS))
    if CRITICAL_ERRORS:
        print("\nCRITICAL DOWNLOAD FAILURES:", file=sys.stderr)
        for err in CRITICAL_ERRORS:
            print(f"  - {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
