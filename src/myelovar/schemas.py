"""Pydantic schemas: pipeline results, funnel, manifest, gene cards.

This module is UI-free by design (pydantic only) so both the API layer and
any frontend consume the same serializable types.
"""
from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Optional

from pydantic import BaseModel, Field


def _utcnow() -> datetime:
    return datetime.now(timezone.utc)


class FunnelStep(BaseModel):
    """Row counts after one pipeline step, with removal reasons."""

    step: str
    label: str
    rows_in: int
    rows_out: int
    removed: int
    removed_reasons: dict[str, int] = Field(default_factory=dict)
    note: str = ""


class RankedVariant(BaseModel):
    """One ranked variant with its full evidence chain (all values optional
    when a data source is missing; explanations omit missing clauses)."""

    variant_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    rank: int = 0
    score: float = 0.0

    # s1 annotation
    consequence_class: str = "noncoding"
    is_splice: bool = False
    af: Optional[float] = None
    dp: Optional[int] = None
    gq: Optional[int] = None
    nearest_gene: Optional[str] = None

    # s2 regulatory
    regulatory_context_score: float = 0.0
    active_tracks: list[str] = Field(default_factory=list)
    in_encode_ccre: bool = False
    in_super_enhancer: bool = False

    # s3 motif
    top_motif_tf: Optional[str] = None
    motif_ref_rel: Optional[float] = None
    motif_alt_rel: Optional[float] = None
    motif_delta: Optional[float] = None
    motif_effect: str = "none"
    motif_disruption_score: float = 0.0

    # s4 gene link
    target_genes: list[str] = Field(default_factory=list)
    link_method: str = "none"
    link_evidence_count: int = 0
    link_confidence: str = "low"
    link_score: float = 0.0
    link_details: list[dict[str, Any]] = Field(default_factory=list)

    # s5 myeloma relevance
    depmap_mm_mean_effect: Optional[float] = None
    depmap_mm_selective: bool = False
    is_known_mm_gene: bool = False
    expressed_in_mm: bool = False
    in_mm_gwas_locus: bool = False
    myeloma_relevance_score: float = 0.0
    gwas_flag: bool = False

    # s6
    explanation: str = ""
    warnings: list[str] = Field(default_factory=list)

    model_config = {"extra": "allow"}


class RunMetadata(BaseModel):
    run_id: str
    mode: str
    created_at: datetime = Field(default_factory=_utcnow)
    config: dict[str, Any] = Field(default_factory=dict)
    manifest_hash: str = ""
    input_description: str = ""
    elapsed_seconds: float = 0.0
    random_seed: int = 0
    warnings: list[str] = Field(default_factory=list)


class SensitivityReport(BaseModel):
    """Weight-perturbation stability of the top-50 ranking."""

    n_perturbations: int
    top_n: int
    jaccard_mean: float
    jaccard_min: float
    spearman_mean: Optional[float] = None
    weight_scheme: str = "each weight perturbed +/-50%, renormalised"


class PipelineResult(BaseModel):
    """High-level return value of run_pipeline(); fully JSON-serializable."""

    schema_version: str = "1.0"
    run: RunMetadata
    funnel: list[FunnelStep] = Field(default_factory=list)
    variants: list[RankedVariant] = Field(default_factory=list)
    total_ranked: int = 0
    sensitivity: Optional[SensitivityReport] = None
    warnings: list[str] = Field(default_factory=list)
    outputs: dict[str, str] = Field(default_factory=dict)


class ManifestEntry(BaseModel):
    """Provenance of one downloaded file (rule: every file traces to a URL)."""

    id: str
    category: str
    source_url: str
    accession: str = ""
    download_date: str
    genome_build: str = "GRCh38"
    cell_type: str = ""
    assay: str = ""
    path: str
    size_bytes: int
    md5: str
    verification: str = ""
    status: str = "ok"          # ok | failed | skipped
    error: str = ""
    notes: str = ""


class DataManifest(BaseModel):
    generated_at: str
    genome_build: str = "GRCh38"
    manifest_hash: str = ""
    resources: list[ManifestEntry] = Field(default_factory=list)


class GeneCard(BaseModel):
    gene: str
    depmap_mm_mean_effect: Optional[float] = None
    depmap_other_mean_effect: Optional[float] = None
    depmap_mm_selective: bool = False
    is_known_mm_gene: bool = False
    known_gene_source: str = ""
    expressed_in_mm: bool = False
    linked_variant_count: int = 0
    top_variants: list[RankedVariant] = Field(default_factory=list)
