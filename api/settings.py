"""Environment-driven API settings (.env.example documents every variable)."""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class Settings:
    data_dir: str = field(default_factory=lambda: os.environ.get("MYELOVAR_DATA_DIR", "data"))
    outputs_dir: str = field(default_factory=lambda: os.environ.get("MYELOVAR_OUTPUTS_DIR", "outputs"))
    demo_dir: str = field(default_factory=lambda: os.environ.get("MYELOVAR_DEMO_DIR", "demo_bundle"))
    validation_dir: str = field(default_factory=lambda: os.environ.get("MYELOVAR_VALIDATION_DIR", "validation"))
    web_dist: str = field(default_factory=lambda: os.environ.get("MYELOVAR_WEB_DIST", "web/dist"))
    max_upload_mb: int = field(default_factory=lambda: int(os.environ.get("MYELOVAR_MAX_UPLOAD_MB", "50")))
    cors_origins: str = field(default_factory=lambda: os.environ.get("MYELOVAR_CORS_ORIGINS", "*"))
    store_patient_data: bool = field(default_factory=lambda: os.environ.get("MYELOVAR_STORE_UPLOADS", "false").lower() == "true")

    @property
    def cors_origin_list(self) -> list[str]:
        return [o.strip() for o in self.cors_origins.split(",") if o.strip()]


def get_settings() -> Settings:
    return Settings()
