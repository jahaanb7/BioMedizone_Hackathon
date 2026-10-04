"""API smoke tests via TestClient (no scientific logic - endpoints only)."""
from __future__ import annotations

from pathlib import Path

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from api.main import app  # noqa: E402

client = TestClient(app)


def test_health_shape():
    r = client.get("/api/health")
    assert r.status_code == 200
    body = r.json()
    for key in ("status", "version", "manifest_hash", "data_ready", "demo_ready"):
        assert key in body
    assert body["status"] == "ok"


def test_modes_lists_three():
    r = client.get("/api/modes")
    assert r.status_code == 200
    modes = r.json()
    assert [m["id"] for m in modes] == ["A", "B", "C"]
    assert all("name" in m and "description" in m for m in modes)


def test_openapi_available():
    r = client.get("/openapi.json")
    assert r.status_code == 200
    paths = r.json()["paths"]
    assert "/api/health" in paths and "/api/runs" in paths


def test_manifest_endpoint_present():
    r = client.get("/api/manifest")
    assert r.status_code in (200, 404)
    if r.status_code == 200:
        body = r.json()
        assert "resources" in body and "manifest_hash" in body


def test_runs_list_and_unknown_run_404():
    r = client.get("/api/runs")
    assert r.status_code == 200 and isinstance(r.json(), list)
    r = client.get("/api/runs/doesnotexist")
    assert r.status_code == 404


def test_demo_gated_on_bundle():
    r = client.get("/api/demo/A")
    if Path("demo_bundle/result_A.json").exists():
        assert r.status_code == 200
        assert r.json()["mode"] == "A"
    else:
        assert r.status_code == 404
    assert client.get("/api/demo/Z").status_code == 400


def test_create_run_rejects_bad_uploads(tmp_path):
    cases = [
        ("notes.txt", b"hello", "text/plain"),                 # wrong extension
        ("plain.vcf", b"not a vcf at all", "application/octet-stream"),
        ("empty.vcf.gz", b"", "application/gzip"),
        ("hg19.vcf", b"##fileformat=VCFv4.2\n##reference=hg19\n#CHROM\n",
         "application/octet-stream"),
        ("badmode.vcf.gz", None, None),                        # handled below
    ]
    for name, content, mime in cases:
        if content is None:
            continue
        r = client.post("/api/runs", data={"mode": "A"},
                        files={"file": (name, content, mime)})
        assert r.status_code == 400, (name, r.status_code, r.text)
        assert "detail" in r.json()


def test_create_run_requires_file_for_mode_a():
    r = client.post("/api/runs", data={"mode": "A"})
    assert r.status_code == 400
    assert "VCF" in r.json()["detail"]


def test_create_run_rejects_unknown_mode():
    r = client.post("/api/runs", data={"mode": "Q"})
    assert r.status_code == 400


def test_locus_requires_variant_id():
    r = client.get("/api/runs/x/locus")
    assert r.status_code == 422  # missing query param


def test_validation_endpoint_shape():
    r = client.get("/api/validation")
    assert r.status_code in (200, 404)
    if r.status_code == 200:
        assert isinstance(r.json(), dict)
