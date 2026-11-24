"""The development app: /v1 (contract C1) and the upstream-compatible legacy routes, over one pool.

    uv run uvicorn irv_metrics.adapters.fastapi:app --port 8001

  /v1/...                          irv_metrics.api: the API the frontend is being ported to
  /metrics/..., /boundaries/...,   irv_metrics.adapters.legacy: what web/ speaks today
  /downloads/...
  /health

Environment: DATABASE_URL; CORS_ORIGINS (default *); IRV_METRICS_DOWNLOADS=0 turns
/v1 downloads off; IRV_METRICS_COUNTRIES=KEN serves one country;
IRV_METRICS_STATEMENT_TIMEOUT_MS (default 30000).
"""
from __future__ import annotations

from contextlib import asynccontextmanager
import os

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from psycopg_pool import ConnectionPool

from ..api import Service, create_router
from ..loaders.common import DEFAULT_DSN, SCHEMA_NAME
from . import legacy

STATEMENT_TIMEOUT_MS = int(os.environ.get("IRV_METRICS_STATEMENT_TIMEOUT_MS", "30000"))

pool = ConnectionPool(os.environ.get("DATABASE_URL", DEFAULT_DSN), open=False, min_size=1, max_size=8, timeout=10,
                      kwargs={"options": f"-c search_path={SCHEMA_NAME},public -c statement_timeout={STATEMENT_TIMEOUT_MS}"})

service = Service(pool, downloads=os.environ.get("IRV_METRICS_DOWNLOADS", "1") != "0",
                  countries={c for c in os.environ.get("IRV_METRICS_COUNTRIES", "").split(",") if c} or None)
legacy_router, catalog = legacy.create_router(pool)


@asynccontextmanager
async def lifespan(_app):
    pool.open(wait=True)
    catalog.refresh()
    service.refresh()
    yield
    pool.close()


app = FastAPI(title="GRI risk data API", version="0.3.0", lifespan=lifespan)
# The web app calls the API directly from the browser (NEXT_PUBLIC_NATIONAL_API_URL).
app.add_middleware(CORSMiddleware, allow_origins=os.environ.get("CORS_ORIGINS", "*").split(","), allow_methods=["GET"],
                   allow_headers=["*"], expose_headers=["ETag"])


@app.get("/health")
def health():
    return {"status": "ok", "series": len(catalog.series), "areas": len(catalog.areas),
            "version": service.state.version if service.state else None}


app.include_router(create_router(service))
app.include_router(legacy_router)
