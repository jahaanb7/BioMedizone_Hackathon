"""Disk-cached HTTP GET/POST JSON client with polite rate limiting."""
from __future__ import annotations

import hashlib
import json
import logging
import threading
import time
from pathlib import Path
from typing import Any, Optional

import requests

log = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path("data/cache/http")
_UA = {"User-Agent": "myelovar/0.1 (research pipeline)"}
_lock = threading.Lock()
_last_call = 0.0


def _throttle(min_interval: float = 0.08) -> None:
    global _last_call
    with _lock:
        wait = min_interval - (time.time() - _last_call)
        if wait > 0:
            time.sleep(wait)
        _last_call = time.time()


def _key(url: str, payload: Optional[dict] = None) -> str:
    h = hashlib.sha1(url.encode())
    if payload:
        h.update(json.dumps(payload, sort_keys=True).encode())
    return h.hexdigest()


def get_json(url: str, cache_dir: Path | str = DEFAULT_CACHE_DIR, *,
             refresh: bool = False, timeout: int = 60, retries: int = 3,
             headers: Optional[dict] = None) -> Any:
    """GET a JSON document, cached on disk (cache_dir/<sha1>.json)."""
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f"{_key(url)}.json"
    if path.exists() and not refresh:
        return json.loads(path.read_text())
    last: Exception | None = None
    for attempt in range(retries):
        _throttle()
        try:
            r = requests.get(url, headers={"Accept": "application/json",
                                           **_UA, **(headers or {})}, timeout=timeout)
            r.raise_for_status()
            data = r.json()
            path.write_text(json.dumps(data))
            return data
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"GET {url} failed after {retries} attempts: {last}")


def post_json(url: str, payload: dict, cache_dir: Path | str = DEFAULT_CACHE_DIR, *,
              refresh: bool = False, timeout: int = 60, retries: int = 3) -> Any:
    """POST JSON, cached by (url, payload) hash."""
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f"{_key(url, payload)}.json"
    if path.exists() and not refresh:
        return json.loads(path.read_text())
    last: Exception | None = None
    for attempt in range(retries):
        _throttle()
        try:
            r = requests.post(url, json=payload,
                              headers={"Accept": "application/json", **_UA},
                              timeout=timeout)
            r.raise_for_status()
            data = r.json()
            path.write_text(json.dumps(data))
            return data
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"POST {url} failed after {retries} attempts: {last}")
