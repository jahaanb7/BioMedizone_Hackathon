"""All REST endpoints. Every operation delegates to the myelovar core or to
disk artifacts produced by it; no scientific logic lives here."""
from __future__ import annotations

import hashlib
import json
import logging
import re
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from fastapi import APIRouter, BackgroundTasks, File, Form, HTTPException, Query, UploadFile
from fastapi.responses import FileResponse, JSONResponse

from api.jobs import STORE
from api.models import (ErrorResponse, GeneCardResponse, HealthResponse,
                        JobStatus, LocusResponse, LocusTrack, ModeInfo,
                        RunCreateResponse, VariantDetail, VariantsPage)
from api.settings import Settings, get_settings
from myelovar import __version__, load_config
from myelovar.reference import ReferenceBundle

log = logging.getLogger("myelovar.api")
router = APIRouter(prefix="/api")

_table_cache: dict[str, tuple[float, pd.DataFrame]] = {}
_bundle_cache: Optional[ReferenceBundle] = None


def _bundle() -> ReferenceBundle:
    global _bundle_cache
    if _bundle_cache is None:
        _bundle_cache = ReferenceBundle(load_config())
    return _bundle_cache


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


# ------------------------------------------------------------- health ----
@router.get("/health", response_model=HealthResponse)
def health() -> HealthResponse:
    s = get_settings()
    try:
        b = _bundle()
        manifest_hash = b.manifest_hash
        data_ready = (Path("data/tracks/tracks.json").exists()
                      and Path("data/data_manifest.json").exists())
    except Exception:  # noqa: BLE001
        manifest_hash, data_ready = "", False
    demo_ready = Path(s.demo_dir, "result_A.json").exists()
    return HealthResponse(status="ok", version=__version__,
                          manifest_hash=manifest_hash,
                          data_ready=data_ready, demo_ready=demo_ready)


@router.get("/manifest")
def manifest() -> dict:
    p = Path(get_settings().data_dir) / "data_manifest.json"
    if not p.exists():
        raise HTTPException(404, "manifest not found - run `make data`")
    return json.loads(p.read_text())


# total_ranked per demo mode, cached by result-file mtime so /modes never
# re-parses a multi-MB result JSON until the bundle actually changes.
_DEMO_TOTAL_CACHE: dict[str, tuple[float, int]] = {}


def _demo_total_ranked(mode: str, demo_dir: str) -> Optional[int]:
    path = Path(demo_dir) / f"result_{mode}.json"
    try:
        mtime = path.stat().st_mtime
    except OSError:
        return None
    hit = _DEMO_TOTAL_CACHE.get(mode)
    if hit is not None and hit[0] == mtime:
        return hit[1]
    try:
        total = int(json.loads(path.read_text()).get("total_ranked") or 0)
    except (OSError, ValueError, json.JSONDecodeError):
        return None
    _DEMO_TOTAL_CACHE[mode] = (mtime, total)
    return total


@router.get("/modes", response_model=list[ModeInfo])
def modes() -> list[ModeInfo]:
    s = get_settings()
    demo_a = Path(s.demo_dir, "result_A.json").exists()
    demo_b = Path(s.demo_dir, "result_B.json").exists()
    total_a = _demo_total_ranked("A", s.demo_dir) if demo_a else None
    total_b = _demo_total_ranked("B", s.demo_dir) if demo_b else None
    return [
        ModeInfo(id="A", name="Whole-genome funnel",
                 description="Real public whole-genome VCF (GIAB HG001/NA12878, GRCh38) "
                             "filtered millions of variants down to a ranked regulatory "
                             "shortlist. Healthy germline genome used to demonstrate the "
                             "pipeline; no real somatic drivers expected.",
                 input_kind="vcf", demo_available=demo_a,
                 demo_run_id="demo-A" if demo_a else None,
                 demo_total_ranked=total_a),
        ModeInfo(id="B", name="Myeloma GWAS locus fine-mapping",
                 description="Candidate causal regulatory variants at published myeloma "
                             "GWAS loci (GWAS Catalog lead SNPs + 1000G EUR LD proxies).",
                 input_kind="loci", demo_available=demo_b,
                 demo_run_id="demo-B" if demo_b else None,
                 demo_total_ranked=total_b),
        ModeInfo(id="C", name="DepMap myeloma somatic mutations (stretch)",
                 description="Somatic mutations in myeloma cell lines from DepMap "
                             "(mostly exome-derived; few non-coding variants expected).",
                 input_kind="depmap", demo_available=False),
    ]


# --------------------------------------------------------------- demo ----
def _demo_payload(mode: str, s: Settings) -> dict:
    mode = mode.upper()
    if mode not in ("A", "B"):
        raise HTTPException(400, f"unknown demo mode {mode}")
    result_path = Path(s.demo_dir) / f"result_{mode}.json"
    if not result_path.exists():
        raise HTTPException(404, f"demo bundle for mode {mode} not built - run `make run` / "
                                 "scripts/make_demo_bundle.py")
    run_id = f"demo-{mode}"
    result = json.loads(result_path.read_text())
    if STORE.get(run_id) is None:
        STORE.create(run_id, mode, result.get("run", {}).get("input_description", "demo"))
        STORE.update(run_id, status="done", fraction=1.0, message="precomputed demo",
                     warnings=result.get("warnings", []),
                     total_ranked=result.get("total_ranked", 0),
                     funnel=result.get("funnel", []),
                     result_available=True, finished_at=_now(),
                     base_dir=str(Path(s.demo_dir) / f"mode_{mode.lower()}"), cache_key=None)
    return {"run_id": run_id, "mode": mode, "result": result}


@router.get("/demo/{mode}")
def demo(mode: str) -> dict:
    return _demo_payload(mode, get_settings())


# ------------------------------------------------------------ runs -------
def _save_upload(upload: UploadFile, dest: Path) -> tuple[int, str]:
    """Stream-upload with size + type + genome-build validation (4xx on bad)."""
    name = (upload.filename or "").lower()
    if not name.endswith((".vcf", ".vcf.gz", ".vcf.bgz")):
        raise HTTPException(400, {"error": f"unsupported file type: {upload.filename!r}",
                                  "hint": "expected .vcf or .vcf.gz"})
    s = get_settings()
    limit = s.max_upload_mb * 1024 * 1024
    dest.parent.mkdir(parents=True, exist_ok=True)
    head = b""
    total = 0
    with open(dest, "wb") as fh:
        while True:
            chunk = upload.file.read(1 << 20)
            if not chunk:
                break
            total += len(chunk)
            if total > limit:
                fh.close()
                dest.unlink(missing_ok=True)
                raise HTTPException(413, f"file exceeds {s.max_upload_mb} MB limit")
            if len(head) < 65536:
                head += chunk[:65536 - len(head)]
            fh.write(chunk)
    if total == 0:
        dest.unlink(missing_ok=True)
        raise HTTPException(400, "empty file")
    if b"##fileformat=VCF" not in head[:200]:
        dest.unlink(missing_ok=True)
        raise HTTPException(400, {"error": "not a VCF file (missing ##fileformat=VCF header)",
                                  "hint": "export a bgzipped or plain VCF on GRCh38/hg38"})
    text_head = head.decode("utf-8", errors="replace")
    if re.search(r"(hg19|GRCh37|b37|hs37d5)", text_head):
        dest.unlink(missing_ok=True)
        raise HTTPException(400, {"error": "input appears to be GRCh37/hg19 - this build is not supported",
                                  "hint": "lift your VCF to GRCh38/hg38 first"})
    return total, hashlib.sha256(head).hexdigest()[:8]


def _execute_run(run_id: str, input_path: Optional[str], mode: str,
                 cache_key: Optional[str]) -> None:
    """Background task: run the core pipeline, updating the job store."""
    from myelovar.pipeline import run_pipeline
    try:
        STORE.update(run_id, status="running")
        cfg = load_config()
        base_out = Path("outputs/runs") / run_id
        cfg.paths.outputs_dir = str(base_out)

        def cb(p: dict) -> None:
            STORE.update(run_id, step=str(p.get("step", "")),
                         fraction=float(p.get("fraction", 0.0)),
                         message=str(p.get("message", "")))

        result = run_pipeline(input=input_path, mode=mode, config=cfg,
                              progress_callback=cb, run_id=run_id)
        base = base_out if mode.upper() == "A" else base_out / f"mode_{mode.lower()}"
        STORE.update(run_id, status="done", fraction=1.0, message="done",
                     warnings=result.warnings, total_ranked=result.total_ranked,
                     funnel=[s.model_dump() for s in result.funnel],
                     result_available=True, finished_at=_now(),
                     base_dir=str(base), cache_key=cache_key)
        log.info("run %s finished: %d ranked", run_id, result.total_ranked)
    except Exception as exc:  # noqa: BLE001 - surfaced to the client verbatim
        log.exception("run %s failed", exc)
        STORE.update(run_id, status="error", error=f"{type(exc).__name__}: {exc}",
                     finished_at=_now())
    finally:
        s = get_settings()
        if input_path and not s.store_patient_data:
            Path(input_path).unlink(missing_ok=True)


@router.post("/runs", response_model=RunCreateResponse)
def create_run(background: BackgroundTasks,
               file: Optional[UploadFile] = File(None),
               mode: str = Form("A")) -> RunCreateResponse:
    mode = mode.upper()
    if mode not in ("A", "B", "C"):
        raise HTTPException(400, f"unknown mode {mode}")
    run_id = uuid.uuid4().hex[:12]

    if mode == "A":
        if file is None:
            raise HTTPException(400, "mode A requires a VCF upload (field 'file')")
        tmp = Path("outputs/uploads") / f"{run_id}.vcf.gz"
        total, _ = _save_upload(file, tmp)
        # run cache: sha256(content) + mode + manifest hash
        digest = hashlib.sha256()
        with open(tmp, "rb") as fh:
            for block in iter(lambda: fh.read(1 << 22), b""):
                digest.update(block)
        try:
            mh = _bundle().manifest_hash
        except Exception:  # noqa: BLE001
            mh = ""
        key = hashlib.sha256(f"{digest.hexdigest()}|{mode}|{mh}".encode()).hexdigest()[:16]
        cached_dir = Path("outputs/runs/byhash") / key
        if (cached_dir / "result.json").exists() or any(
                cached_dir.glob("**/result.json")):
            STORE.create(run_id, mode, f"uploaded VCF ({total} bytes, cached)")
            base = next(cached_dir.glob("**/result.json")).parent
            result = json.loads((base / "result.json").read_text())
            STORE.update(run_id, status="done", fraction=1.0,
                         message="served from cache", warnings=result.get("warnings", []),
                         total_ranked=result.get("total_ranked", 0),
                         funnel=result.get("funnel", []),
                         result_available=True, finished_at=_now(),
                         base_dir=str(base), cache_key=key)
            tmp.unlink(missing_ok=True)
            return RunCreateResponse(run_id=run_id, status="done", mode=mode,
                                     input_description="uploaded VCF (cached result)",
                                     cached=True)
        # keep a content-addressed copy for cache reuse
        keep = Path("outputs/runs/byhash") / key
        keep.mkdir(parents=True, exist_ok=True)
        final_in = keep / "input.vcf.gz"
        tmp.rename(final_in)
        input_path = str(final_in)
        desc = f"uploaded VCF ({total} bytes)"
    elif mode == "B":
        input_path = None
        desc = "Mode B: GWAS-locus candidates (loci from cache)"
        key = None
    else:
        input_path = None
        desc = "Mode C: DepMap somatic mutations in myeloma lines"
        key = None

    STORE.create(run_id, mode, desc)
    background.add_task(_execute_run, run_id, input_path, mode, key)
    return RunCreateResponse(run_id=run_id, status="queued", mode=mode,
                             input_description=desc, cached=False)


def _synth_job_from_disk(run_id: str) -> Optional[JobStatus]:
    """Rebuild job status from outputs after an API restart."""
    base = Path("outputs/runs") / run_id
    candidates = [base, base / "mode_b", base / "mode_c"]
    for cand in candidates:
        if (cand / "result.json").exists():
            result = json.loads((cand / "result.json").read_text())
            mode = result.get("run", {}).get("mode", "A")
            job = STORE.create(run_id, mode, result.get("run", {}).get("input_description", ""))
            STORE.update(job.run_id, status="done", fraction=1.0, message="restored from disk",
                         warnings=result.get("warnings", []),
                         total_ranked=result.get("total_ranked", 0),
                         funnel=result.get("funnel", []),
                         result_available=True, finished_at=_now(),
                         base_dir=str(cand))
            return STORE.get(run_id)
    return None


def _get_job(run_id: str) -> JobStatus:
    job = STORE.get(run_id) or _synth_job_from_disk(run_id)
    if job is None:
        raise HTTPException(404, f"unknown run_id {run_id}")
    return job


def _ranked_csv(job: JobStatus) -> Path:
    if not job.base_dir:
        raise HTTPException(409, "run has not finished yet")
    p = Path(job.base_dir) / "ranked_variants.csv"
    if not p.exists():
        raise HTTPException(404, "ranked_variants.csv missing for this run")
    return p


def _load_table(job: JobStatus) -> pd.DataFrame:
    csv = _ranked_csv(job)
    key = str(csv)
    mtime = csv.stat().st_mtime
    hit = _table_cache.get(key)
    if hit and hit[0] == mtime:
        return hit[1]
    df = pd.read_csv(csv, low_memory=False)
    _table_cache[key] = (mtime, df)
    return df


@router.get("/runs", response_model=list[JobStatus])
def list_runs() -> list[JobStatus]:
    return STORE.all()


@router.get("/runs/{run_id}", response_model=JobStatus)
def run_status(run_id: str) -> JobStatus:
    return _get_job(run_id)


@router.get("/runs/{run_id}/funnel")
def run_funnel(run_id: str) -> dict:
    job = _get_job(run_id)
    return {"run_id": run_id, "mode": job.mode, "steps": job.funnel,
            "status": job.status}


def _filter_table(df: pd.DataFrame, *, min_score: Optional[float] = None,
                  gene: Optional[str] = None, tf: Optional[str] = None,
                  chrom: Optional[str] = None, confidence: Optional[str] = None,
                  in_super_enhancer: Optional[bool] = None,
                  locus_id: Optional[str] = None, motif_effect: Optional[str] = None) -> pd.DataFrame:
    out = df
    if min_score is not None:
        out = out[out["score"] >= min_score]
    if gene:
        pat = rf"(^|,)\s*{re.escape(gene)}\s*(,|$)"
        out = out[out["target_genes"].fillna("").str.contains(pat, regex=True)
                  | out["nearest_gene"].fillna("").eq(gene)]
    if tf:
        out = out[out["top_motif_tf"].fillna("").str.upper() == tf.upper()]
    if chrom:
        c = chrom if chrom.startswith("chr") else f"chr{chrom}"
        out = out[out["chrom"] == c]
    if confidence:
        out = out[out["link_confidence"] == confidence]
    if in_super_enhancer is not None:
        col = out["in_super_enhancer"].fillna(False).astype(bool)
        out = out[col == in_super_enhancer]
    if locus_id and "locus_id" in out.columns:
        out = out[out["locus_id"] == locus_id]
    if motif_effect and "motif_effect" in out.columns:
        out = out[out["motif_effect"] == motif_effect]
    return out


@router.get("/runs/{run_id}/variants", response_model=VariantsPage)
def run_variants(run_id: str, min_score: Optional[float] = None,
                 gene: Optional[str] = None, tf: Optional[str] = None,
                 chrom: Optional[str] = None, confidence: Optional[str] = None,
                 in_super_enhancer: Optional[bool] = None,
                 locus_id: Optional[str] = None, motif_effect: Optional[str] = None,
                 sort: str = Query("score", pattern="^(score|rank|pos|variant_id)$"),
                 order: str = Query("desc", pattern="^(asc|desc)$"),
                 limit: int = Query(50, ge=1, le=500), offset: int = Query(0, ge=0)) -> VariantsPage:
    job = _get_job(run_id)
    df = _load_table(job)
    df = _filter_table(df, min_score=min_score, gene=gene, tf=tf, chrom=chrom,
                       confidence=confidence, in_super_enhancer=in_super_enhancer,
                       locus_id=locus_id, motif_effect=motif_effect)
    if sort == "pos":
        df = df.assign(_k=df["pos"])
    else:
        df = df.assign(_k=df[sort] if sort in df.columns else df["score"])
    df = df.sort_values("_k", ascending=(order == "asc"), kind="stable")
    total = len(df)
    page = df.iloc[offset:offset + limit].drop(columns=["_k"])
    rows = json.loads(page.to_json(orient="records"))
    return VariantsPage(run_id=run_id, total=total, offset=offset, limit=limit, rows=rows)


def _find_variant(df: pd.DataFrame, variant_id: str) -> dict:
    hit = df[df["variant_id"] == variant_id]
    if hit.empty:
        raise HTTPException(404, f"variant {variant_id} not in this run's ranked table")
    rec = json.loads(hit.head(1).to_json(orient="records"))[0]
    for k, v in list(rec.items()):
        if isinstance(v, str) and k in ("target_genes", "active_tracks"):
            rec[k] = [x for x in v.split(",") if x]
        if isinstance(v, str) and k == "warnings":
            rec[k] = [x for x in v.split("|") if x]
        if isinstance(v, str) and k == "link_details" and v:
            try:
                rec[k] = json.loads(v)
            except json.JSONDecodeError:
                rec[k] = []
    return rec


@router.get("/runs/{run_id}/variants/{variant_id}", response_model=VariantDetail)
def run_variant(run_id: str, variant_id: str) -> VariantDetail:
    job = _get_job(run_id)
    df = _load_table(job)
    rec = _find_variant(df, variant_id)
    return VariantDetail(run_id=run_id, variant=rec,
                         explanation=rec.get("explanation", ""),
                         warnings=rec.get("warnings", []))


def _demo_locus_path(variant_id: str) -> Optional[Path]:
    safe = re.sub(r"[^A-Za-z0-9]+", "_", variant_id)
    p = Path(get_settings().demo_dir) / "loci" / f"{safe}.json"
    return p if p.exists() else None


@router.get("/runs/{run_id}/locus", response_model=LocusResponse)
def run_locus(run_id: str, variant_id: str = Query(...),
              window: int = Query(50_000, ge=1000, le=500_000)) -> LocusResponse:
    job = _get_job(run_id)
    # demo runs serve precomputed locus payloads (fully offline)
    if run_id.startswith("demo-"):
        p = _demo_locus_path(variant_id)
        if p:
            data = json.loads(p.read_text())
            data["cached"] = True
            return LocusResponse(**data)
        raise HTTPException(404, f"no precomputed locus payload for {variant_id} in the demo bundle")

    df = _load_table(job)
    rec = _find_variant(df, variant_id)
    try:
        return _build_locus_payload(rec, variant_id, window)
    except FileNotFoundError as exc:
        raise HTTPException(503, f"reference data unavailable for locus view: {exc}")


def _build_locus_payload(rec: dict, variant_id: str, window: int) -> LocusResponse:
    """Build a locus payload from a normalized variant record. Shared by the
    live endpoint and scripts/make_demo_bundle.py (precomputed demo loci)."""
    chrom, pos = rec["chrom"], int(rec["pos"])
    start, end = max(0, pos - window), pos + window
    b = _bundle()
    try:
        tracks_out: list[LocusTrack] = []
        for t in b.peak_tracks():
            idx = b.track_index(t["name"])
            arr = idx.chroms.get(chrom)
            if arr is None:
                continue
            starts, ends = arr
            lo = int(np.searchsorted(starts, start, side="left"))
            hits = []
            i = lo
            while i < len(starts) and starts[i] < end and len(hits) < 200:
                if ends[i] > start:
                    hits.append([int(starts[i]), int(ends[i])])
                i += 1
            if hits:
                tracks_out.append(LocusTrack(
                    name=t["name"], display_name=t.get("display_name", ""),
                    cell_type=t.get("cell_type", ""), assay=t.get("assay", ""),
                    category=t.get("category", ""), intervals=hits))
        ccre_idx = b.ccre()
        ccre_count = len(ccre_idx.overlapping_intervals(chrom, start, end))
        genes_df = b.genes()
        sub = genes_df[genes_df["chrom"] == chrom]
        genes_out = []
        if not sub.empty:
            s_arr = sub["start"].to_numpy()
            e_arr = sub["end"].to_numpy()
            sel = np.where((e_arr > start) & (s_arr < end))[0][:200]
            for i in sel:
                row = sub.iloc[int(i)]
                genes_out.append({"gene": row["gene"], "start": int(row["start"]),
                                  "end": int(row["end"]), "strand": row["strand"],
                                  "tss": int(row["tss"]), "gene_type": row["gene_type"]})
        links = rec.get("link_details") or []
        if isinstance(links, str):
            links = json.loads(links) if links else []
        motif = None
        if rec.get("top_motif_tf"):
            try:
                m = b.motifs()[rec["top_motif_tf"]]
                W = b.config.motif.window
                fa = b.fasta(chrom)
                ref_seq = fa.fetch(chrom, max(0, pos - 1 - W), pos - 1 + len(rec["ref"]) + W).upper()
                alt_seq = ref_seq[:W] + rec["alt"].upper() + ref_seq[W + len(rec["ref"]):]
                motif = {"tf": rec["top_motif_tf"], "matrix_id": m["matrix_id"],
                         "effect": rec.get("motif_effect", "none"),
                         "ref_rel": rec.get("motif_ref_rel"), "alt_rel": rec.get("motif_alt_rel"),
                         "delta": rec.get("motif_delta"),
                         "ref_seq": ref_seq, "alt_seq": alt_seq, "variant_offset": W,
                         "pwm": m["pwm"].tolist()}
            except FileNotFoundError:
                motif = {"tf": rec["top_motif_tf"], "error": "reference FASTA unavailable"}
        return LocusResponse(variant_id=variant_id, chrom=chrom, pos=pos,
                             window=[start, end], tracks=tracks_out,
                             ccre_count=ccre_count, genes=genes_out,
                             links=links, motif=motif, cached=False)
    except FileNotFoundError:
        raise


@router.get("/runs/{run_id}/export.csv")
def run_export(run_id: str) -> FileResponse:
    job = _get_job(run_id)
    csv = _ranked_csv(job)
    return FileResponse(csv, media_type="text/csv",
                        filename=f"myelovar_{run_id}_ranked_variants.csv")


@router.get("/genes/{symbol}", response_model=GeneCardResponse)
def gene_card(symbol: str) -> GeneCardResponse:
    b = _bundle()
    try:
        depmap = b.depmap()
    except FileNotFoundError:
        depmap = pd.DataFrame(columns=["gene"])
    row = depmap[depmap["gene"].str.upper() == symbol.upper()]
    mm = b.mm_genes()
    known = mm[mm["gene"].str.upper() == symbol.upper()] if len(mm) else mm
    top: list[dict] = []
    linked = 0
    for job in STORE.all():
        if not job.result_available:
            continue
        try:
            df = _load_table(job)
        except HTTPException:
            continue
        pat = rf"(^|,)\s*{re.escape(symbol)}\s*(,|$)"
        hit = df[df["target_genes"].fillna("").str.contains(pat, regex=True)]
        linked += len(hit)
        top.extend(json.loads(hit.nlargest(5, "score").to_json(orient="records")))
    top.sort(key=lambda r: -r.get("score", 0))
    return GeneCardResponse(
        gene=symbol.upper(),
        depmap_mm_mean_effect=float(row["depmap_mm_mean_effect"].iloc[0]) if len(row) else None,
        depmap_other_mean_effect=float(row["depmap_other_mean_effect"].iloc[0]) if len(row) else None,
        depmap_mm_selective=bool(row["depmap_mm_selective"].iloc[0]) if len(row) else False,
        is_known_mm_gene=bool(len(known)),
        known_gene_source=str(known["source"].iloc[0]) if len(known) else "",
        expressed_in_mm=bool(row["expressed_in_mm"].iloc[0]) if len(row) else False,
        linked_variant_count=linked, top_variants=top[:10])


@router.get("/validation")
def validation() -> JSONResponse:
    p = Path(get_settings().validation_dir) / "results.json"
    if not p.exists():
        raise HTTPException(404, "validation results not built - run scripts/run_validation.py")
    return JSONResponse(json.loads(p.read_text()))


@router.get("/validation/figures/{name}")
def validation_figure(name: str) -> FileResponse:
    if "/" in name or ".." in name:
        raise HTTPException(400, "bad figure name")
    p = Path(get_settings().validation_dir) / "figures" / name
    if not p.exists():
        raise HTTPException(404, f"figure {name} not found")
    return FileResponse(p)
