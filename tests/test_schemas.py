"""JSON round-trips for the pydantic schemas + funnel/rank invariants of real outputs."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from myelovar.schemas import (DataManifest, FunnelStep, ManifestEntry, PipelineResult,
                              RankedVariant, RunMetadata, SensitivityReport)
from myelovar import load_result


def _sample_result() -> PipelineResult:
    return PipelineResult(
        run=RunMetadata(run_id="t1", mode="A", manifest_hash="abc123",
                        input_description="test input"),
        funnel=[FunnelStep(step="s1", label="annotate", rows_in=100, rows_out=60,
                           removed=40, removed_reasons={"common_af": 40})],
        variants=[RankedVariant(variant_id="chr1:100:A:T", chrom="chr1", pos=100,
                                ref="A", alt="T", rank=1, score=0.9)],
        total_ranked=1,
        sensitivity=SensitivityReport(n_perturbations=5, top_n=10,
                                      jaccard_mean=0.8, jaccard_min=0.5),
        warnings=["w1"],
        outputs={"csv": "outputs/ranked_variants.csv"},
    )


def test_pipeline_result_json_round_trip():
    res = _sample_result()
    again = PipelineResult.model_validate_json(res.model_dump_json())
    assert again == res
    # plain-dict path too
    assert PipelineResult.model_validate(res.model_dump()) == res


def test_load_result_from_disk(tmp_path):
    p = tmp_path / "result.json"
    res = _sample_result()
    p.write_text(res.model_dump_json(indent=2))
    loaded = load_result(p)
    assert loaded == res
    with pytest.raises(FileNotFoundError):
        load_result(tmp_path / "missing.json")


def test_manifest_round_trip_and_hash_field():
    m = DataManifest(
        generated_at="2026-10-04T00:00:00+00:00",
        manifest_hash="deadbeef",
        resources=[ManifestEntry(id="giab", category="variants",
                                 source_url="https://example.org/x.vcf.gz",
                                 download_date="2026-10-04", path="data/x.vcf.gz",
                                 size_bytes=10, md5="0" * 32, status="ok")],
    )
    again = DataManifest.model_validate_json(m.model_dump_json())
    assert again == m
    assert again.resources[0].id == "giab"


def test_funnel_steps_serialize_with_reasons():
    s = FunnelStep(step="s1", label="x", rows_in=10, rows_out=5, removed=5,
                   removed_reasons={"a": 3, "b": 2}, note="n")
    d = json.loads(s.model_dump_json())
    assert d["removed_reasons"] == {"a": 3, "b": 2}
    assert d["removed"] == d["rows_in"] - d["rows_out"]


# ------------------------------------------- invariants of real outputs ----
def _real_result_path() -> Path | None:
    for cand in (Path("outputs/result.json"), Path("outputs/mode_b/result.json"),
                 Path("demo_bundle/result_A.json"), Path("demo_bundle/result_B.json")):
        if cand.exists():
            return cand
    return None


def test_real_output_funnel_and_ranking_invariants():
    """Runs whenever a pipeline output exists (make run / demo build)."""
    p = _real_result_path()
    if p is None:
        pytest.skip("no pipeline output yet - run `make run` first")
    res = load_result(p)
    assert res.funnel, "funnel must not be empty"
    # strict funnel invariant: rows never grow after load, removed == difference
    for step in res.funnel:
        assert step.rows_out <= step.rows_in, f"{step.step} grew rows"
        assert step.removed == step.rows_in - step.rows_out, \
            f"{step.step} removed={step.removed} != rows_in-rows_out"
    # rows are monotonically non-increasing down the funnel
    for prev, step in zip(res.funnel, res.funnel[1:]):
        assert step.rows_in == prev.rows_out, "funnel stages must chain"
    # ranked variants: contiguous ranks, non-increasing scores
    if res.variants:
        assert [v.rank for v in res.variants] == list(range(1, len(res.variants) + 1))
        scores = [v.score for v in res.variants]
        assert all(a >= b - 1e-9 for a, b in zip(scores, scores[1:]))
    # JSON round-trip of the whole thing
    assert PipelineResult.model_validate_json(res.model_dump_json()) == res
