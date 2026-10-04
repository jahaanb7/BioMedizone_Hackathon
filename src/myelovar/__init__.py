"""Myelovar: myeloma-specific regulatory variant prioritization (research only).

Core library with zero UI/web dependencies. Typical use:

    from myelovar import run_pipeline, load_config
    result = run_pipeline("input.vcf.gz", mode="A", config=load_config())
"""
from myelovar.config import AppConfig, load_config
from myelovar.pipeline import load_result, run_pipeline
from myelovar.schemas import PipelineResult

__version__ = "0.1.0"

__all__ = ["run_pipeline", "load_result", "load_config", "AppConfig",
           "PipelineResult", "__version__"]
