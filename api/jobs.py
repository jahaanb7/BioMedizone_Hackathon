"""In-process job store with a swappable interface.

The interface (submit/ get/ update) is deliberately tiny so Celery/RQ can
replace it later without touching routes or science code.
"""
from __future__ import annotations

import threading
from datetime import datetime, timezone
from typing import Optional

from api.models import JobStatus


class JobStore:
    """Thread-safe registry of run jobs (dict-backed; swap for Redis/Celery)."""

    def __init__(self) -> None:
        self._jobs: dict[str, JobStatus] = {}
        self._lock = threading.Lock()

    def create(self, run_id: str, mode: str, input_description: str) -> JobStatus:
        job = JobStatus(run_id=run_id, mode=mode, status="queued",
                        input_description=input_description,
                        created_at=datetime.now(timezone.utc).isoformat())
        with self._lock:
            self._jobs[run_id] = job
        return job

    def get(self, run_id: str) -> Optional[JobStatus]:
        with self._lock:
            return self._jobs.get(run_id)

    def update(self, run_id: str, **fields) -> JobStatus:
        with self._lock:
            job = self._jobs.get(run_id)
            if job is None:
                raise KeyError(f"unknown run_id {run_id}")
            for k, v in fields.items():
                setattr(job, k, v)
            return job

    def all(self) -> list[JobStatus]:
        with self._lock:
            return list(self._jobs.values())


STORE = JobStore()
