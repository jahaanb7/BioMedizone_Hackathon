"""Pydantic request/response models for the REST API (validation only)."""
from __future__ import annotations

from typing import Any, Optional

from pydantic import BaseModel, Field


class HealthResponse(BaseModel):
    status: str
    version: str
    manifest_hash: str
    data_ready: bool
    demo_ready: bool


class ModeInfo(BaseModel):
    id: str
    name: str
    description: str
    input_kind: str
    demo_available: bool
    demo_run_id: Optional[str] = None
    # ranked-variant count of the precomputed demo result (0/None when absent);
    # lets the UI show real counts on a cold start without loading the full
    # multi-MB demo payload.
    demo_total_ranked: Optional[int] = None


class JobStatus(BaseModel):
    run_id: str
    mode: str
    status: str                       # queued | running | done | error
    step: str = ""
    fraction: float = 0.0
    message: str = ""
    input_description: str = ""
    created_at: str = ""
    finished_at: Optional[str] = None
    error: Optional[str] = None
    warnings: list[str] = Field(default_factory=list)
    total_ranked: int = 0
    funnel: list[dict[str, Any]] = Field(default_factory=list)
    result_available: bool = False
    base_dir: Optional[str] = None
    cache_key: Optional[str] = None


class RunCreateResponse(BaseModel):
    run_id: str
    status: str
    mode: str
    input_description: str
    cached: bool = False


class VariantsPage(BaseModel):
    run_id: str
    total: int
    offset: int
    limit: int
    rows: list[dict[str, Any]]


class VariantDetail(BaseModel):
    run_id: str
    variant: dict[str, Any]
    explanation: str
    warnings: list[str]


class LocusTrack(BaseModel):
    name: str
    display_name: str = ""
    cell_type: str = ""
    assay: str = ""
    category: str = ""
    intervals: list[list[int]] = Field(default_factory=list)


class LocusResponse(BaseModel):
    variant_id: str
    chrom: str
    pos: int
    window: list[int]
    tracks: list[LocusTrack] = Field(default_factory=list)
    ccre_count: int = 0
    genes: list[dict[str, Any]] = Field(default_factory=list)
    links: list[dict[str, Any]] = Field(default_factory=list)
    motif: Optional[dict[str, Any]] = None
    cached: bool = False


class GeneCardResponse(BaseModel):
    gene: str
    depmap_mm_mean_effect: Optional[float] = None
    depmap_other_mean_effect: Optional[float] = None
    depmap_mm_selective: bool = False
    is_known_mm_gene: bool = False
    known_gene_source: str = ""
    expressed_in_mm: bool = False
    linked_variant_count: int = 0
    top_variants: list[dict[str, Any]] = Field(default_factory=list)


class ErrorResponse(BaseModel):
    detail: str
    hint: Optional[str] = None
