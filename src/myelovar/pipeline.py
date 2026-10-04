"""run_pipeline(): step orchestration, progress reporting, output writing.

The pipeline is decoupled from any UI: progress goes through an optional
progress_callback, logs go through logging, results are typed pydantic objects
that serialize cleanly to JSON.
"""
from __future__ import annotations

import json
import logging
import time
import uuid
from pathlib import Path
from typing import Any, Callable, Optional

import pandas as pd

from myelovar.config import AppConfig, load_config
from myelovar.reference import ReferenceBundle
from myelovar.schemas import (FunnelStep, PipelineResult, RankedVariant,
                              RunMetadata, SensitivityReport)
from myelovar.steps.s0_load_input import step_load_input
from myelovar.steps.s1_annotate_filter import step_annotate_filter
from myelovar.steps.s2_regulatory_overlap import step_regulatory_overlap
from myelovar.steps.s3_motif_effect import step_motif_effect
from myelovar.steps.s4_gene_link import step_gene_link
from myelovar.steps.s5_myeloma_relevance import step_myeloma_relevance
from myelovar.steps.s6_score_rank import step_score_rank

log = logging.getLogger(__name__)

ProgressFn = Callable[[dict[str, Any]], None]

# pipeline step boundaries -> progress fraction (for UI progress bars)
_PROGRESS = {"s0": 0.05, "s1": 0.20, "s2": 0.40, "s3": 0.65,
             "s4": 0.80, "s5": 0.88, "s6": 0.97, "done": 1.0}


def _emit(cb: Optional[ProgressFn], step: str, message: str) -> None:
    if cb is not None:
        cb({"step": step, "fraction": _PROGRESS.get(step, 0.0),
            "message": message})


def _outputs_base(config: AppConfig) -> Path:
    base = Path(config.paths.outputs_dir)
    mode = config.mode.upper()
    if mode != "A":
        base = base / f"mode_{mode.lower()}"
    base.mkdir(parents=True, exist_ok=True)
    return base


def _records_to_variants(df: pd.DataFrame, cap: int) -> list[RankedVariant]:
    """Convert the top-`cap` ranked rows to JSON-clean RankedVariant models."""
    sub = df.head(cap)
    raw = json.loads(sub.to_json(orient="records"))
    out: list[RankedVariant] = []
    for rec in raw:
        tg = rec.get("target_genes") or ""
        rec["target_genes"] = [g for g in str(tg).split(",") if g]
        w = rec.get("warnings") or ""
        rec["warnings"] = [x for x in str(w).split("|") if x]
        at = rec.get("active_tracks") or ""
        rec["active_tracks"] = [x for x in str(at).split(",") if x]
        ld = rec.get("link_details")
        if isinstance(ld, str) and ld:
            try:
                rec["link_details"] = json.loads(ld)
            except json.JSONDecodeError:
                rec["link_details"] = []
        for k in ("in_super_enhancer", "in_encode_ccre", "gwas_flag",
                  "depmap_mm_selective", "is_known_mm_gene", "expressed_in_mm",
                  "in_mm_gwas_locus", "is_splice"):
            if k in rec and rec[k] is not None:
                rec[k] = bool(rec[k])
        out.append(RankedVariant.model_validate(rec))
    return out


def run_pipeline(input: Optional[str] = None, mode: str = "A",
                 config: Optional[AppConfig] = None,
                 progress_callback: Optional[ProgressFn] = None,
                 *, max_records: Optional[int] = None,
                 run_id: Optional[str] = None) -> PipelineResult:
    """Run the full prioritisation pipeline; returns a serializable result.

    Args:
        input: path to a VCF (Mode A); ignored in Mode B (loci built from cache).
        mode: "A" (whole-genome funnel), "B" (GWAS loci) or "C" (DepMap somatic).
        config: AppConfig; loaded from config.yaml when omitted.
        progress_callback: called with {"step", "fraction", "message"} dicts.
        max_records: optional cap on input records (fast demos/tests).
        run_id: optional caller-supplied id (API generates one).
    """
    cfg = config or load_config()
    cfg = cfg.model_copy(deep=True)
    cfg.mode = mode.upper()
    run_id = run_id or uuid.uuid4().hex[:12]
    t0 = time.time()
    base = _outputs_base(cfg)
    step_dir = base / "step_outputs"
    step_dir.mkdir(parents=True, exist_ok=True)
    bundle = ReferenceBundle(cfg)

    funnel: list[FunnelStep] = []
    warnings: list[str] = []

    _emit(progress_callback, "s0", "Loading input variants")
    df, f0 = step_load_input(cfg, input_path=input, max_records=max_records)
    total_input = len(df)
    funnel.append(f0)

    _emit(progress_callback, "s1", "Quality / blacklist / common-variant / coding filters")
    df, f1, w = step_annotate_filter(df, bundle=bundle, config=cfg,
                                     out_dir=step_dir, mode=cfg.mode)
    funnel.append(f1); warnings += w

    _emit(progress_callback, "s2", "Overlapping active regulatory elements")
    df, f2, w = step_regulatory_overlap(df, bundle=bundle, config=cfg, out_dir=step_dir)
    funnel.append(f2); warnings += w

    _emit(progress_callback, "s3", "Scanning TF motif disruption")
    df, f3, w = step_motif_effect(df, bundle=bundle, config=cfg, out_dir=step_dir)
    funnel.append(f3); warnings += w

    _emit(progress_callback, "s4", "Linking enhancers to target genes")
    df, f4, w = step_gene_link(df, bundle=bundle, config=cfg, out_dir=step_dir)
    funnel.append(f4); warnings += w

    _emit(progress_callback, "s5", "Scoring myeloma relevance")
    df, f5, w = step_myeloma_relevance(df, bundle=bundle, config=cfg, out_dir=step_dir)
    funnel.append(f5); warnings += w

    _emit(progress_callback, "s6", "Ranking and writing explanations")
    df, f6, sens, w = step_score_rank(df, bundle=bundle, config=cfg, out_dir=step_dir,
                                      total_input_rows=total_input)
    funnel.append(f6); warnings += w

    # ------------------------------------------------------------ outputs --
    funnel_payload = {
        "run_id": run_id, "mode": cfg.mode, "input_total": total_input,
        "final_ranked": len(df),
        "steps": [s.model_dump() for s in funnel],
    }
    (base / "funnel.json").write_text(json.dumps(funnel_payload, indent=2))
    if sens is not None:
        (base / "sensitivity.json").write_text(json.dumps(sens.model_dump(), indent=2))

    cap = cfg.runtime.max_explained_rows
    variants = _records_to_variants(df, cap)
    meta = RunMetadata(
        run_id=run_id, mode=cfg.mode,
        config={
            "filter": cfg.filter.model_dump(),
            "motif": cfg.motif.model_dump(),
            "link": cfg.link.model_dump(),
            "score_weights": cfg.score_weights,
            "regulatory_context_weights": cfg.regulatory_context_weights,
            "myeloma_relevance_weights": cfg.myeloma_relevance_weights,
            "random_seed": cfg.runtime.random_seed,
        },
        manifest_hash=bundle.manifest_hash,
        input_description=f0.note,
        elapsed_seconds=round(time.time() - t0, 3),
        random_seed=cfg.runtime.random_seed,
        warnings=warnings,
    )
    result = PipelineResult(
        run=meta, funnel=funnel, variants=variants, total_ranked=len(df),
        sensitivity=sens, warnings=warnings,
        outputs={
            "ranked_csv": str(base / "ranked_variants.csv"),
            "funnel_json": str(base / "funnel.json"),
            "result_json": str(base / "result.json"),
        },
    )
    (base / "result.json").write_text(result.model_dump_json(indent=2))
    _emit(progress_callback, "done",
          f"Done: {total_input} -> {len(df)} ranked variants "
          f"({result.run.elapsed_seconds}s)")
    log.info("run %s complete: %s -> %d ranked in %.1fs",
             run_id, cfg.mode, len(df), result.run.elapsed_seconds)
    return result


def load_result(path: str | Path = "outputs/result.json") -> PipelineResult:
    """Load a previously written result.json."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"result file not found: {p}")
    return PipelineResult.model_validate_json(p.read_text())
