# GRI Infra-Risk-Vis API

REST API, written in Python, includes database definition and Tileserver.

The features API is a light wrapper over tables stored in Postgres, using the
PostGIS extension for spatial data types.

The tiles API calls into Terracotta to read cloud-optimised geotiffs and serve
PNG tiles.

Tiles must be loaded separately - there are no endpoints for ingesting data.

## Installation / Build

### Docker

```bash
docker-compose -f docker-compose-prod-build.yaml build backend
```

### Running Locally

```bash
# cd to this backend app directory
cd backend

# create a virtual environment (using venv or another method if you prefer)
python -m venv venv
source venv/bin/activate
# irv-metrics (the national metrics API) first: the backend depends on it by
# name and it is not on PyPI
pip install -e ../irv-metrics
pip install -e .[dev]

# run the application
uvicorn backend.app.main:app --port 8888 --reload
```

## Database migrations

The schema of the backend database is defined by the SQLAlchemy models in
`backend/db/models.py` and versioned with [Alembic](https://alembic.sqlalchemy.org/):

- `alembic.ini` and `migrations/` live next to `pyproject.toml` and are copied
  into the Docker image at `/code/alembic.ini` and `/code/migrations`.
- `migrations/env.py` connects with the same `PG*` environment variables as
  the app (no URL in `alembic.ini`).
- `migrations/versions/` holds one file per revision.
  - `0001_baseline` creates the tables in `models.py` and runs `CREATE EXTENSION
    IF NOT EXISTS postgis`.
  - `0002_metrics` sets up the `irv-metrics` schema
- Autogenerate only compares the default (`public`) schema, and ignores tables
  owned by database extensions (PostGIS's `spatial_ref_sys`, pg_cron's tables)
  and the `topology`, `tiger`, `tiger_data`, `cron` and `metrics` schemas.

Run these from this directory, with `PG*` set for the target database:

```bash
alembic upgrade head      # apply all migrations (what the `migrate` service runs)
alembic current           # show the database's revision
alembic check             # fail if models.py and the database differ
alembic upgrade head --sql  # print the SQL instead of running it
```

### Existing databases

For already-existing databases, run `alembic stamp 0001_baseline` once, then use
`upgrade head` as usual.

### Changing the schema

For each model change, make a new alembic revision:

1. Edit `backend/db/models.py`.
2. Bring a dev database up to date (`alembic upgrade head`), then run
   `alembic revision --autogenerate -m "short description"`.
3. Read the generated file in `migrations/versions/` and fix it by hand where
   autogenerate falls short (renames might show up as drop + add; data
   migrations may be needed, server-side defaults and custom SQL have to be
   written). Geometry columns render as `geoalchemy2.Geometry(...)`; decide
   whether to allow alembic to create the spatial index or to be more explicit.
4. `alembic upgrade head`, then `alembic check` should report no changes, and
   `alembic downgrade -1` should work.
5. Commit the revision with the model change.

Schema objects that are not in `models.py` (e.g. a SQL-defined schema) go in
hand-written revisions (`alembic revision -m "..."`, then `op.execute(...)`).
If they live in their own schema, add it to `EXCLUDED_SCHEMAS` in
`migrations/env.py` so autogenerate leaves it alone

For example, `0002_metrics` is hand-written: it creates the `metrics` schema and
runs the installed `irv-metrics` package's `irv_metrics/schema/*.sql` (via
`importlib.resources`).

Downgrading to `0001_baseline` drops the schema and its data. The metric
package's loaders leave the schema to Alembic. Any change to those SQL files
needs a new revision..

### Reset

`alembic downgrade base && alembic upgrade head` drops every table in the
migrations and all their data, then recreates them empty. The `dev-reset-db`
compose service does this for the dev stack. Don't run it on production.

### Tests

The tests check the migrations against a real database: `upgrade head` on an
empty database, then autogenerate must find no differences with `models.py`;
`downgrade base` then `upgrade head` must work, via the API and the CLI. They
need a **throwaway** PostGIS database (its GRI tables are dropped) given as
`DATABASE_URL`, and are skipped without it:

```bash
pip install -e ../irv-metrics
pip install -e .[dev]
docker run -d --name gri-test-db -p 5499:5432 \
    -e POSTGRES_DB=gri_test -e POSTGRES_USER=gri_test -e POSTGRES_PASS=password \
    kartoza/postgis:14-3.1
# wait until it accepts connections (initialisation restarts postgres, ~10 s)
until docker logs gri-test-db 2>&1 | grep -q "restarting in foreground"; do sleep 1; done; sleep 3
DATABASE_URL=postgresql+psycopg2://gri_test:password@localhost:5499/gri_test pytest
docker rm -f gri-test-db
```

`tests/test_metrics_migration.py` checks `0002_metrics` (needs irv-metrics
installed). CI runs these and irv-metrics' tests
(`.github/workflows/backend-tests.yml`).

## Metrics API

This mounts the irv-metrics `/v1` router at `/metrics`.

Environment variables:
- `IRV_METRICS_DATABASE_URL`
- `IRV_METRICS_POOL_MIN_SIZE` defaults to 1
- `IRV_METRICS_POOL_MAX_SIZE` defaults to 4
- `IRV_METRICS_POOL_TIMEOUT_S` defaults to 10
- `IRV_METRICS_STATEMENT_TIMEOUT_MS` defaults to 30000
- `IRV_METRICS_DOWNLOADS` defaults to 1; 0 turns CSV downloads off
- `IRV_METRICS_COUNTRIES` defaults to every loaded country

See [../irv-metrics/README.md](../irv-metrics/README.md) for the API and
loaders.

### API performance

The `benchmark_pixel_driller.py` script runs queries against `localhost` running
the API to roughly benchmark performance:

```bash
$ python benchmark_pixel_driller.py -c 4 -n 100
Requests: 100 | Concurrency: 4
Mean latency: 1543.54 ms
Total wall time (approx): 154.35 s
P95 latency: 2276.02 ms
P99 latency: 2416.50 ms
Status counts:
  200: 100

$  python benchmark_pixel_driller.py -c 8 -n 100
Requests: 100 | Concurrency: 8
Mean latency: 3179.33 ms
Total wall time (approx): 317.93 s
P95 latency: 4366.23 ms
P99 latency: 4833.38 ms
Status counts:
  200: 100
```

## Configuration

Environment variables:

- use `.env` to define environment variables
- use [`PG*`](https://www.postgresql.org/docs/current/libpq-envars.html) to
  define database connection details. See `.env.example` for an example

```
PYTHONPATH=/code/backend
# Features API
PGDATABASE=
PGUSER=
PGPASSWORD=
PGHOST=

# Tiles API
RASTER_BASE_PATH=/data # Path at-which raster tiles can be found (must match the Tiles-db loaded path)
TC_DRIVER_PATH=postgresql://foo:bar@tiles-db # Tiles-db MySQL Host (__NOTE__: Does not require database in the URL - this is parsed internally.)
TC_ALLOWED_ORIGINS_METADATA='["*"]'
TC_ALLOWED_ORIGINS_TILES='["*"]'
TC_PNG_COMPRESS_LEVEL=0
TC_RESAMPLING_METHOD="nearest"
TC_REPROJECTION_METHOD="nearest"
```

### Tileserver

Tileserver endpoints `/tiles/*` require Cloud-Optimised GeoTIFFs to be mounted
at the same path as they have been ingested into the database. The tileserver
endpoints require raster data (tiffs) to have been pre-loaded using the `etl`
pipelines.

The base path for rasters must be set using the environment variable
`RASTER_BASE_PATH`

**NOTE**: TC_DRIVER_PATH is not used internally - for Terracotta the path is
built programatically based on the URL

Tileserver also provides information about each tile database. This metadata is
stored in the `raster_tile_sources` table in the main `global_prod`/`global_dev`
database.

#### Colormaps

Colormaps for use with Terracotta can be generated using the `/colormap`
endpoint.

#### Categorical Data

Categorical rasters are supported. Categorical colormaps can either be included
in the `config.py`, or passed with each tile request, as per the Terracotta
documentation:
https://terracotta-python.readthedocs.io/en/latest/tutorials/categorical.html

**NOTE**: Only `{pixel :(RGBA)}` explicit color maps are supported in either the
request or `config,py`

If included in `config.py` the key must match an existing MySQL database, with
pre-loaded categorical raster(s).

e.g. for a `land_cover` database the entries would be similar to the following:

```json
CATEGORICAL_COLOR_MAPS = {
    "land_cover": {
        0: (0, 0, 0, 255),
        10: (255, 255, 100, 255),
        11: (255, 255, 100, 255),
    }
}
```
