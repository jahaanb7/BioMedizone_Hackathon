"""s6: weighted score, ranking, explanations and weight sensitivity analysis."""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from myelovar.config import AppConfig
from myelovar.explain import build_explanation, row_warnings
from myelovar.reference import ReferenceBundle
from myelovar.schemas import FunnelStep, SensitivityReport

log = logging.getLogger(__name__)

COMPONENTS = ("regulatory_context", "motif_disruption", "link",
              "myeloma_relevance", "gwas")
COMPONENT_COLS = ("regulatory_context_score", "motif_disruption_score",
                  "link_score", "myeloma_relevance_score", "_gwas_comp")


def _component_matrix(df: pd.DataFrame) -> np.ndarray:
    return np.column_stack([
        df["regulatory_context_score"].to_numpy(float),
        df["motif_disruption_score"].to_numpy(float),
        df["link_score"].to_numpy(float),
        df["myeloma_relevance_score"].to_numpy(float),
        df["gwas_flag"].astype(float).to_numpy(),
    ])


def sensitivity_analysis(X: np.ndarray, base_scores: np.ndarray, *,
                         n_perturbations: int, seed: int,
                         top_n: int = 50) -> SensitivityReport:
    """Perturb each weight +/-50% (renormalised); report top-N stability."""
    rng = np.random.default_rng(seed)
    base_top = set(np.argsort(-base_scores)[:top_n].tolist())
    n = len(base_scores)
    sub = rng.choice(n, size=min(50_000, n), replace=False) \
        if n > 60_000 else np.arange(n)
    base_rank_sub = np.argsort(np.argsort(-base_scores[sub]))
    jaccards: list[float] = []
    spearmans: list[float] = []
    for _ in range(n_perturbations):
        pert = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        pert *= rng.uniform(0.5, 1.5, size=5)
        pert = pert / pert.sum()
        scores = X @ pert
        top = set(np.argsort(-scores)[:top_n].tolist())
        inter = len(top & base_top)
        union = len(top | base_top) or 1
        jaccards.append(inter / union)
        pert_rank_sub = np.argsort(np.argsort(-scores[sub]))
        if base_rank_sub.std() > 0 and pert_rank_sub.std() > 0:
            c = np.corrcoef(base_rank_sub, pert_rank_sub)[0, 1]
            spearmans.append(float(c))
    return SensitivityReport(
        n_perturbations=n_perturbations, top_n=top_n,
        jaccard_mean=float(np.mean(jaccards)), jaccard_min=float(np.min(jaccards)),
        spearman_mean=float(np.mean(spearmans)) if spearmans else None,
        weight_scheme="each weight perturbed +/-50% then renormalised; "
                      f"Spearman computed on a {len(sub)}-row subsample")


def step_score_rank(df: pd.DataFrame, *, bundle: ReferenceBundle,
                    config: AppConfig, out_dir: Path,
                    total_input_rows: int) -> tuple[pd.DataFrame, FunnelStep,
                                                    SensitivityReport | None, list[str]]:
    """Score, rank, explain; writes outputs/ranked_variants.csv (+ parquet)."""
    rows_in = len(df)
    warnings: list[str] = []
    df = df.copy()

    weights = config.normalised_score_weights()
    # X columns follow COMPONENTS order (config keys are identical)
    w_vec = np.array([weights[k] for k in COMPONENTS])
    X = _component_matrix(df)
    df["score"] = np.round(X @ w_vec, 6)

    df = df.sort_values(["score", "variant_id"], ascending=[False, True]).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)

    # sensitivity analysis (seeded; row-aligned with the sorted frame)
    sens: SensitivityReport | None = None
    if len(df) >= 10:
        sens = sensitivity_analysis(_component_matrix(df),
                                    df["score"].to_numpy(float),
                                    n_perturbations=config.runtime.sensitivity_perturbations,
                                    seed=config.runtime.random_seed)

    # explanations + per-row warnings (real values only)
    total = len(df)
    expl: list[str] = []
    warn_col: list[list[str]] = []
    for rec in df.to_dict(orient="records"):
        expl.append(build_explanation(rec, int(rec["rank"]), total))
        warn_col.append(row_warnings(rec))
    df["explanation"] = expl
    df["warnings"] = ["|".join(w) for w in warn_col]

    csv_path = Path(config.paths.outputs_dir)
    if config.mode.upper() == "B":
        csv_path = csv_path / "mode_b"
    csv_path.mkdir(parents=True, exist_ok=True)
    export = df.copy()
    export.to_csv(csv_path / "ranked_variants.csv", index=False)
    df.to_parquet(out_dir / "s6_ranked.parquet", index=False)

    step = FunnelStep(
        step="s6", label="Score, rank, explain", rows_in=rows_in, rows_out=len(df),
        removed=0,
        note=f"weights={{{', '.join(f'{k}: {v:.2f}' for k, v in weights.items())}}}; "
             + (f"top-50 Jaccard mean={sens.jaccard_mean:.2f}" if sens else ""))
    log.info("s6: ranked %d variants; top-50 sensitivity Jaccard=%.3f",
             len(df), sens.jaccard_mean if sens else float("nan"))
    return df, step, sens, warnings
