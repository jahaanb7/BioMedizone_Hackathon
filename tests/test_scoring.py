"""s6 scoring: weight normalisation, deterministic ranking, sensitivity, explanations."""
from __future__ import annotations

import numpy as np
import pandas as pd

from myelovar.config import AppConfig
from myelovar.steps.s6_score_rank import COMPONENTS, sensitivity_analysis, step_score_rank


def _make_frame(n: int = 12, seed: int = 5) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    df = pd.DataFrame({
        "variant_id": [f"chr1:{1000 + i}:A:T" for i in range(n)],
        "chrom": ["chr1"] * n,
        "pos": [1000 + i for i in range(n)],
        "ref": ["A"] * n,
        "alt": ["T"] * n,
        "regulatory_context_score": rng.uniform(0, 1, n),
        "motif_disruption_score": rng.uniform(0, 1, n),
        "link_score": rng.uniform(0, 1, n),
        "myeloma_relevance_score": rng.uniform(0, 1, n),
        "gwas_flag": rng.uniform(0, 1, n) > 0.7,
        # columns the explainer reads (kept realistic but sparse)
        "nearest_gene": ["KRAS"] * n,
        "target_genes": ["KRAS,NRAS"] * n,
        "link_method": ["ABC"] * n,
        "link_confidence": ["medium"] * n,
        "active_tracks": ["H3K27ac_ENCSR758OEC"] * n,
        "motif_effect": ["broken"] * n,
        "top_motif_tf": ["MYC"] * n,
        "motif_ref_rel": [0.9] * n,
        "motif_alt_rel": [0.5] * n,
    })
    return df.reset_index(drop=True)


def _run(df, tmp_path, tag):
    cfg = AppConfig(mode="A")
    cfg.paths.outputs_dir = str(tmp_path / f"out_{tag}")
    out = tmp_path / f"work_{tag}"
    out.mkdir()
    return step_score_rank(df, bundle=None, config=cfg, out_dir=out,
                           total_input_rows=len(df))


def test_normalised_weights_sum_to_one():
    cfg = AppConfig()
    w = cfg.normalised_score_weights()
    assert abs(sum(w.values()) - 1.0) < 1e-12
    assert set(w) == set(COMPONENTS)


def test_scores_are_weighted_sum_and_ranks_are_sorted(tmp_path):
    df = _make_frame()
    cfg = AppConfig(mode="A")
    cfg.paths.outputs_dir = str(tmp_path / "out")
    out = tmp_path / "work"
    out.mkdir()
    ranked, step, sens, warnings = step_score_rank(
        df, bundle=None, config=cfg, out_dir=out, total_input_rows=len(df))

    w = cfg.normalised_score_weights()
    cols = ["regulatory_context_score", "motif_disruption_score", "link_score",
            "myeloma_relevance_score", "gwas_flag"]
    expected = sum(ranked[c].to_numpy(float) * w[k]
                   for c, k in zip(cols, COMPONENTS))
    # the step rounds scores to 6 decimals
    assert np.allclose(ranked["score"].to_numpy(), np.round(expected, 6))
    assert list(ranked["rank"]) == list(range(1, len(ranked) + 1))
    scores = ranked["score"].to_numpy()
    assert np.all(np.diff(scores) <= 1e-9)          # non-increasing
    assert step.rows_out == step.rows_in == len(df) and step.removed == 0


def test_scoring_is_deterministic_across_runs(tmp_path):
    df = _make_frame()
    r1, _, s1, _ = _run(df, tmp_path, "a")
    r2, _, s2, _ = _run(df, tmp_path, "b")
    assert r1["score"].tolist() == r2["score"].tolist()
    assert r1["rank"].tolist() == r2["rank"].tolist()
    assert r1["explanation"].tolist() == r2["explanation"].tolist()
    # sensitivity uses a fixed seed -> identical stability estimates
    assert s1 is not None and s2 is not None
    assert s1.jaccard_mean == s2.jaccard_mean
    assert s1.jaccard_min <= s1.jaccard_mean <= 1.0
    assert 0.0 <= s1.jaccard_min


def test_sensitivity_analysis_seed_reproducible():
    rng = np.random.default_rng(3)
    X = rng.uniform(0, 1, size=(80, 5))
    base = X @ np.array([0.25, 0.2, 0.2, 0.25, 0.1])
    a = sensitivity_analysis(X, base, n_perturbations=20, seed=1729, top_n=20)
    b = sensitivity_analysis(X, base, n_perturbations=20, seed=1729, top_n=20)
    assert a.jaccard_mean == b.jaccard_mean and a.jaccard_min == b.jaccard_min
    assert a.spearman_mean is None or -1.0 <= a.spearman_mean <= 1.0


def test_explanations_use_only_real_values(tmp_path):
    df = _make_frame(10)
    ranked, _, _, _ = _run(df, tmp_path, "expl")
    for text in ranked["explanation"]:
        assert isinstance(text, str) and len(text) > 40
        # no template placeholders or invented values may leak through
        assert "{" not in text and "None" not in text and "nan" not in text
        # every explanation anchors on its variant and carries real values
        assert text.startswith("chr") and ("MM.1S" in text or "Rank" in text)
