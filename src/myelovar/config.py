"""Configuration loading. The core never reads globals: every consumer passes
an AppConfig (loaded from config.yaml or constructed from arguments)."""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Optional

import yaml
from pydantic import BaseModel, Field


class PathsConfig(BaseModel):
    data_dir: str = "data"
    outputs_dir: str = "outputs"
    cache_dir: str = "data/cache"
    demo_dir: str = "demo_bundle"
    manifest: str = "data/data_manifest.json"
    track_registry: str = "data/tracks/tracks.json"


class FilterConfig(BaseModel):
    min_dp: int = 10
    min_gq: int = 20
    max_af: float = 0.01
    remove_common: bool = True
    autosomes_only: bool = True


class MotifConfig(BaseModel):
    window: int = 30
    min_strong: float = 0.80
    min_drop: float = 0.20


class LinkConfig(BaseModel):
    distance_kb: int = 100
    decay_scale_kb: int = 10


class DepmapConfig(BaseModel):
    mm_effect_threshold: float = -0.5
    selectivity_margin: float = 0.30
    expressed_pct_threshold: int = 50


class GwasConfig(BaseModel):
    window_kb: int = 50
    ld_r2: float = 0.6
    population: str = "1000GENOMES:phase_3:EUR"
    efo_trait_query: str = "multiple myeloma"


class ApiConfig(BaseModel):
    max_upload_mb: int = 50
    cors_origins: str = "*"
    host: str = "127.0.0.1"
    port: int = 8000


class RuntimeConfig(BaseModel):
    random_seed: int = 1729
    sensitivity_perturbations: int = 60
    max_explained_rows: int = 5000


class AppConfig(BaseModel):
    mode: str = "A"
    input: Optional[str] = None
    paths: PathsConfig = Field(default_factory=PathsConfig)
    filter: FilterConfig = Field(default_factory=FilterConfig)
    motif: MotifConfig = Field(default_factory=MotifConfig)
    link: LinkConfig = Field(default_factory=LinkConfig)
    regulatory_context_weights: dict[str, float] = Field(
        default_factory=lambda: {
            "accessibility": 0.30, "h3k27ac": 0.25, "h3k4me1": 0.10,
            "h3k4me3": 0.10, "super_enhancer": 0.15, "tf_peak": 0.10,
        }
    )
    myeloma_relevance_weights: dict[str, float] = Field(
        default_factory=lambda: {
            "depmap": 0.40, "known_gene": 0.30, "expressed": 0.15, "gwas": 0.15,
        }
    )
    score_weights: dict[str, float] = Field(
        default_factory=lambda: {
            "regulatory_context": 0.25, "motif_disruption": 0.20,
            "link": 0.20, "myeloma_relevance": 0.25, "gwas": 0.10,
        }
    )
    depmap: DepmapConfig = Field(default_factory=DepmapConfig)
    gwas: GwasConfig = Field(default_factory=GwasConfig)
    api: ApiConfig = Field(default_factory=ApiConfig)
    runtime: RuntimeConfig = Field(default_factory=RuntimeConfig)

    def normalised_score_weights(self) -> dict[str, float]:
        """Score weights, renormalised to sum to 1 (guards hand-edited configs)."""
        total = sum(self.score_weights.values()) or 1.0
        return {k: v / total for k, v in self.score_weights.items()}


def load_config(path: Optional[str | Path] = None, overrides: Optional[dict[str, Any]] = None) -> AppConfig:
    """Load YAML config; environment variables MYELOVAR_* override paths."""
    cfg_path = Path(path or os.environ.get("MYELOVAR_CONFIG", "config.yaml"))
    data: dict[str, Any] = {}
    if cfg_path.exists():
        with open(cfg_path, "r", encoding="utf-8") as fh:
            loaded = yaml.safe_load(fh) or {}
        if not isinstance(loaded, dict):
            raise ValueError(f"Config {cfg_path} must be a mapping, got {type(loaded)}")
        data = loaded
    if overrides:
        data = {**data, **{k: v for k, v in overrides.items() if v is not None}}
    cfg = AppConfig.model_validate(data)
    if os.environ.get("MYELOVAR_DATA_DIR"):
        cfg.paths.data_dir = os.environ["MYELOVAR_DATA_DIR"]
    if os.environ.get("MYELOVAR_OUTPUTS_DIR"):
        cfg.paths.outputs_dir = os.environ["MYELOVAR_OUTPUTS_DIR"]
    if os.environ.get("MYELOVAR_MAX_UPLOAD_MB"):
        cfg.api.max_upload_mb = int(os.environ["MYELOVAR_MAX_UPLOAD_MB"])
    if os.environ.get("MYELOVAR_CORS_ORIGINS"):
        cfg.api.cors_origins = os.environ["MYELOVAR_CORS_ORIGINS"]
    return cfg
