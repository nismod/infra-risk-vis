"""
Metrics: the irv-metrics `/v1` API, mounted at `/metrics`.

Behind the proxy this is `/api/metrics/v1/...` (traefik strips `/api`).

The routes, response models, errors and caching headers all come from the
`irv_metrics` package (../irv-metrics); this module only provides the database
pool and the settings.

The data lives in the `metrics` schema, created by the Alembic revision
`0002_metrics` and filled by the `metrics-load` compose service.

The connection pool is psycopg3 (the rest of the backend currently uses
SQLAlchemy over psycopg2) and connects with the same `PG*` variables, unless
IRV_METRICS_DATABASE_URL is set.

Environment variables used here:

- IRV_METRICS_DATABASE_URL: libpq conninfo or URL (default: empty, i.e. PG*)
- IRV_METRICS_POOL_MIN_SIZE (default 1), IRV_METRICS_POOL_MAX_SIZE (default 4):
  per backend process. With SQLAlchemy's default pool (5 + 10 overflow) and
  terracotta, one process stays well below PostgreSQL's default max_connections
  (100).
- IRV_METRICS_POOL_TIMEOUT_S (default 10): wait for a free connection, then 503
  `busy`.
- IRV_METRICS_STATEMENT_TIMEOUT_MS (default 30000): slower queries are cancelled
  and answered 503 `busy`.
- IRV_METRICS_DOWNLOADS (default 1): 0 turns the CSV downloads off (404
  `downloads_disabled`; the catalogue reports `download: false`).
- IRV_METRICS_COUNTRIES (default: all loaded): comma-separated ISO3 codes to
  serve.
"""

from contextlib import asynccontextmanager
import logging
from os import getenv

import irv_metrics.adapters.legacy as legacy
from irv_metrics.api import Service, create_router
from psycopg_pool import ConnectionPool

logger = logging.getLogger(__name__)

SCHEMA = "metrics"
DATABASE_URL = getenv("IRV_METRICS_DATABASE_URL", "")
POOL_MIN_SIZE = int(getenv("IRV_METRICS_POOL_MIN_SIZE", "1"))
POOL_MAX_SIZE = int(getenv("IRV_METRICS_POOL_MAX_SIZE", "4"))
POOL_TIMEOUT_S = float(getenv("IRV_METRICS_POOL_TIMEOUT_S", "10"))
STATEMENT_TIMEOUT_MS = int(getenv("IRV_METRICS_STATEMENT_TIMEOUT_MS", "30000"))
DOWNLOADS = getenv("IRV_METRICS_DOWNLOADS", "1") != "0"
COUNTRIES = {c.strip().upper() for c in getenv("IRV_METRICS_COUNTRIES", "").split(",") if c.strip()} or None

pool = ConnectionPool(
    DATABASE_URL,
    open=False,
    min_size=POOL_MIN_SIZE,
    max_size=POOL_MAX_SIZE,
    timeout=POOL_TIMEOUT_S,
    name="metrics",
    kwargs={
        "application_name": "gri-backend metrics",
        "options": f"-c search_path={SCHEMA},public -c statement_timeout={STATEMENT_TIMEOUT_MS}",
    },
)

service = Service(pool, downloads=DOWNLOADS, countries=COUNTRIES)
router = create_router(service)
legacy_router = legacy.create_router(pool)


@asynccontextmanager
async def lifespan(_app):
    """Open the pool at startup, close it at shutdown.

    Startup does not wait for the database: the service reads its catalogue on
    the first request if the database is not reachable (or not loaded) yet.
    """
    pool.open(wait=False)
    try:
        service.refresh()
    except Exception as error:  # the rest of the API must still start
        logger.warning("national metrics not ready at startup (%s); will retry on first request", error)
    try:
        yield
    finally:
        pool.close()
