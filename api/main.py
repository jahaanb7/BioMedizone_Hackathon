"""Thin FastAPI application: mounts the router, CORS, and the built React SPA.

No scientific logic lives here (or anywhere in api/) - everything delegates to
the importable myelovar core or to disk artifacts it produced.
"""
from __future__ import annotations

import logging
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles

from api.routes import router
from api.settings import get_settings

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(name)s %(message)s")

settings = get_settings()

app = FastAPI(
    title="Myelovar API",
    version="0.1.0",
    description="Myeloma regulatory variant prioritization - thin REST layer "
                "over the myelovar core package.",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.cors_origin_list,
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(router)


@app.get("/api")
def api_root() -> dict:
    return {"docs": "/docs", "openapi": "/openapi.json", "health": "/api/health"}


@app.get("/")
def index():
    # serve the SPA directly (no redirect: a redirect would land on
    # /index.html, which client-side routing must also handle)
    if (_dist / "index.html").exists():
        return FileResponse(_dist / "index.html")
    return RedirectResponse("/docs")


# ----------------------------------------------------- static SPA ------
_dist = Path(settings.web_dist)
if _dist.exists():
    # hashes of built assets (immutable, cache forever); index.html never cached
    app.mount("/assets", StaticFiles(directory=_dist / "assets"), name="assets")

    @app.get("/{path:path}", include_in_schema=False)
    def spa(path: str) -> FileResponse:
        cand = _dist / path
        if path and cand.is_file():
            return FileResponse(cand)
        return FileResponse(_dist / "index.html")
else:

    @app.get("/{path:path}", include_in_schema=False)
    def spa_missing(path: str) -> FileResponse:
        from fastapi import HTTPException
        raise HTTPException(
            404, "web UI not built - run `make web` (or `npm run build` in web/)")
