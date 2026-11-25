# irv-metrics: national metrics data cube and `/v1` API

This is a Python package (`irv_metrics`) for interacting with a PostgreSQL +
PostGIS data cube of risk-related metrics.

It includes loaders that pick up CSVs (driven by a YAML dataset config)
and defines the `/v1` API contract [`openapi.json`](openapi.json).

The backend server application (`../backend`) imports `irv_metrics` and mounts
its router at `/metrics` (`backend/app/routers/metrics.py`).

The `metrics-load` compose service (profile `metrics`) runs the CSV loader with
`--config config/gri`.

## Tests

Install dev dependencies and run tests:

```sh
pip install -e .[dev]
pytest # run unit tests; skip database tests
DATABASE_URL=postgresql://user:pass@localhost:5432/db pytest
```

With `DATABASE_URL` the loader, boundary-set, overlay, standalone and
schema-guard tests run on their fixtures, each in a transaction that is rolled
back.

The `/v1` tests against a loaded cube (`tests/test_api_v1.py`) need KEN loaded
from the scrape, and `tests/test_config_sources.py` needs the published CSVs
(`CSV_DATA`); they skip otherwise. CI (`.github/workflows/backend-tests.yml`)
runs both suites on `kartoza/postgis:14-3.1`.

## Layout

```
irv_metrics/
  specs.py               upstream API (routes, parameters, JSON layouts) <-> canonical cube, one description per dataset
  schema/0001_init.sql   schema definition
  loaders/common.py      DB helpers shared by loaders
  loaders/csv.py         published CSVs + boundary files -> cube, driven by the dataset config
  config.py              dataset config: typed catalogue from config/gri/*.yaml
  queries.py             SQL queries the API needs
  api/models.py          /v1 response models and parameter types
  api/service.py         /v1 reads without HTTP: catalogue, coverage, values, areas, curves, boundaries, CSV
  api/routes.py          /v1 FastAPI router (create_router); errors, ETag/304 and Cache-Control live here
  adapters/legacy.py     upstream-compatible routes used until the /v1 port
  adapters/fastapi.py    dev app: /v1 + legacy routes + /health over one pool
```

Design and evidence: [docs/data-model/README.md](../docs/data-model/README.md).

## Schema

Everything lives in the Postgres schema `metrics` (GRI hosts it under that name).
`schema/0001_init.sql` is applied when `metrics.observation` is missing.

- **Boundary sets.** Every area belongs to a `boundary_set` of kind `admin`,
  `overlay` or `group`. `area` is keyed `(boundary_set, area_id)`, where
  `area_id` is the published id (ISO3 at admin level 0, the CGAZ shapeID below),
  with an integer `area_key` that facts are keyed on.
- **Country scope.** `area_country(area_key, iso3, share)` relates areas to
  countries: one row with share 1 per admin area; one row per overlapped
  country for overlay areas, with the fraction of the area inside it
  (`loaders.common.relate_overlay`, by intersection with level 0 of the set's
  `countries_from` admin set, dropping overlaps below `min_share`); one row per
  member for groups. Country-scoped queries filter through it, so a basin
  shared by Kenya and Uganda is listed, with the same values, under both.
- **Slim facts.** `observation(series_id int4, area_key int4, value)`, primary
  key `(series_id, area_key)` plus an index on `area_key`. No country and no
  load id per row. Not partitioned.
- **Load records.** `load_record` has one row per `(dataset, boundary_set,
  iso3)`, replaced on reload; `iso3` is null for a set loaded globally. The API
  reloads its in-memory catalogue when a load finishes.
- **Structure from the config.** Sections, units, dimensions and members,
  boundary sets and levels, datasets, stored metrics and dataset dimensions are
  seeded from the YAML dataset config (`loaders.common.seed_catalogue`), under
  the config's (C2) ids, by every loader before it writes values; existing rows
  are updated to match. Nothing is deleted: after removing a dataset or member
  from the config, rebuild with `--reset`. Derived `share_of` metrics are not
  stored (the API computes them). A series the config does not allow fails the
  load. Any config directory works, e.g.
  `--config datasets/examples/standalone-jam` (`tests/test_standalone.py`).

## Run

Without Docker, with any PostGIS database:

```sh
uv sync --group dev
export DATABASE_URL=postgresql://gri:gri@localhost:5433/gri
```

Load data:

```sh
uv run python -m irv_metrics.loaders.csv --data-dir ../data [--config config/gri] [--countries KEN MOZ] \
    [--if-loaded skip|replace] [--reset] [--report load.json]
```

Every file pattern in the config (`source.files` of each dataset, `levels[].files`
of each boundary set) is relative to `--data-dir`. On the dev machine the data dir
is a directory of symlinks (not committed; create it once, at the repository root):

```sh
mkdir -p data/boundaries
ln -s /home/mert2014/projects/pawehe/csv-data/KEN data/KEN    # published CSVs, one directory per country
ln -s /home/mert2014/projects/pawehe/csv-data/MOZ data/MOZ
ln -s /data/projects/svi/incoming_data/boundaries data/boundaries/cgaz   # geoBoundariesCGAZ_ADM{0,1,2}.gpkg
```

The CGAZ GeoPackages are global; only the loaded countries' features are read
(with `pyogrio`, which bundles GDAL). Countries and admin levels are discovered
from the file names. As in the original ETL, a blank admin-0 file is skipped as a
stub, and a dataset with admin-1/2 files but no usable admin-0 file fails the
load. Rows published without a value are stored as NULL. Each country is loaded
in one transaction (its areas, then one load record per dataset); a row the
source rules cannot place, or a value on an unknown area, rolls the country back
and the loader exits 1. It prints, per country and dataset, files, rows, rows
placed, skipped (by the rules' reason) and empty, series and time.
`--if-loaded skip` (default) leaves a loaded country alone; `replace` removes
its areas and values and loads it again.

Some stored metrics are computed at load rather than published: `COMPUTED` in
`loaders/csv.py` declares them per dataset id. For GRI, the tropical-cyclone wind
bands (`area_in_band`, `area_in_band_share`) are the differences of consecutive
thresholds (`area_exceeding[_share]`; the top band is the top threshold),
pairing `wind_threshold` and `wind_band` members in order. The load prints them
as `computed`. (C2 has no construct for this yet; the declaration is checked
against the config when the loader starts.)

A single-country config whose files carry no `{iso3}` (the `standalone-jam`
example) is loaded for the countries listed in its vocabulary.

### Comparing two cubes

```sh
uv run python -m irv_metrics.compare postgresql://.../gri_csv postgresql://.../gri_scrape \
    --rounded a:3 --split-dim demographic_group --wind-bands
```

Merge-joins both cubes' series, values (relative tolerance 1e-9) and curve
points, compares areas, metrics and dimension members, and lists differences
grouped by dataset with examples. `--rounded a:3` counts differences no larger
than A's rounding to 3 decimals (the CSVs' precision) apart; `--wind-bands`
checks B's wind bands against differences of A's thresholds. Exit 0 when
identical.

### From the scrape

The scrape follows the upstream API's names (`specs.py`); it writes the config's
ids (the upstream `concentration-curves` path is
`adaptation-outcomes/concentration-curves`) and does not store the curves'
upstream `hazard=flood`. It seeds the structure from `--config` (default
`config/gri`), which must define `cgaz-6.0.0`. A scraped unit that differs from
the config's is a warning (the legacy routes render the config's units); a
series the config does not allow fails the country (exit 1).

The scrape loader prints, per country, how many facts it found, how many were reached
through more than one response, any conflicts between them, any JSON leaves the
mapping could not place, and whether published percentages match the derivation
rule. It stops nothing on a warning; read the output. After a reload the API picks
up the new data within 5 seconds.

## API

### `/v1` (contract C1)

The schema is [`openapi.json`](openapi.json) (frozen; generated from `contract/v1.py`),
the reasoning [docs/api-v1.md](../docs/api-v1.md). Paths are relative to the mount point.

```
GET /v1/countries                                          listCountries
GET /v1/countries/{iso3}                                   getCountry: sets, levels, coverage
GET /v1/datasets                                           getCatalogue (from config/gri, not the database)
GET /v1/countries/{iso3}/boundaries/{level}?set=           getBoundaries (application/geo+json)
GET /v1/countries/{iso3}/datasets/{section}/{dataset}/values?level=&set=&metric=&{dim}=…    getValues
GET /v1/countries/{iso3}/datasets/{section}/{dataset}/areas?area_id=…&set=&metric=&{dim}=…  getAreaObservations
GET /v1/countries/{iso3}/datasets/{section}/{dataset}/curve?area_id=&set=&metric=&{dim}=…   getCurve
GET /v1/countries/{iso3}/datasets/{section}/{dataset}/download?level=&set=                 downloadDataset (CSV)
```

```sh
curl 'localhost:8001/v1/countries/KEN/datasets/risk/capital-stock/areas?area_id=KEN&area_id=32016919B10079274489584&metric=loss&sector=total'
curl 'localhost:8001/v1/countries/KEN/datasets/exposure/population/values?level=1&metric=population_share&demographic_group=female'
```

- **Catalogue and data.** `getCatalogue` is built from the YAML config
  (`irv_metrics.config.load_catalogue`); the database supplies the countries,
  coverage, series ids and values. Dataset ids are the config's: the scrape
  stores the curves as `concentration-curves` with a `hazard` member, which
  `/v1` maps to `adaptation-outcomes/concentration-curves` and ignores.
- **Filters.** Every query parameter other than `level`, `metric`, `set` and
  `area_id` is a dimension. `values` and `curve` need the complete assignment of
  one series (dimensions the metric is not keyed by are ignored); `areas` takes
  a partial one. Derived metrics (`derive.share_of`) are computed on read when
  `metric=` names one; `aggregate_of` metrics are not served yet.
- **Errors** are `{"detail": {"code", "message"}}`: 404 `unknown_country`,
  `unknown_dataset`, `unknown_area`, `unknown_level`, `unknown_boundary_set`,
  `downloads_disabled`; 422 `invalid_filters` (including parameter validation);
  503 `busy` when the pool is exhausted or a statement times out.
- **Caching.** The catalogue `version` is a hash of the catalogue and the latest
  finished load; it is every response's `ETag`, and `If-None-Match` gets a 304.
  `Cache-Control: public, max-age=3600` on countries, catalogue and boundaries,
  `max-age=300` on data, `no-store` on errors. The service notices a finished
  load within 5 seconds.
- **Download CSV**: `country_iso3, boundary_set, level, area_id, area_name,
  [share], {dimensions…}, metric, unit, [x], value`, one row per area and stored
  series (per point for curves), values at full precision.

**Mounting elsewhere (GRI).** The package does not import `contract/`. Everything
host-specific is in the router, so a host only needs a pool:

```python
from irv_metrics.api import Service, create_router
service = Service(pool)                    # psycopg_pool.ConnectionPool with search_path=metrics,public
app.include_router(create_router(service), prefix="/api/national-metrics")
# after pool.open(): service.refresh()
```

`Service(pool, catalogue=None, countries={"KEN"}, downloads=False)` serves one
country or turns downloads off. Set a `statement_timeout` on the pool's
connections so slow queries become 503 `busy`.

### Legacy routes

The upstream-compatible routes (`adapters/legacy.py`) are what the original web
app (`web/`) and its end-to-end test speak until the `/v1` port (F3) deletes them:

```
GET /boundaries/{iso3}/{level}
GET /metrics/{family}/{layer}/{iso3}/{level}?...
GET /metrics/{family}/{layer}/{iso3}/profile/{adm_id}?...
GET /metrics/{family}/{layer}/{iso3}/ladder/{adm_id}?...
GET /metrics/concentration-curves/{iso3}?...
GET /downloads/{ISO3}/{ISO3}_adm{1|2}_{card_slug}_metrics.csv
```

Downloads are generated from the cube (one row per area and series, canonical
dimension columns). The original API served the published research CSVs as-is;
those files are not available here, so the column layout differs from them.
CORS is open by default (`CORS_ORIGINS`) because the web app calls the API from the browser.
Unknown parameters or combinations without data return 422 `{"detail": {"code": "invalid_filters", ...}}`.

Dev app settings: `IRV_METRICS_DOWNLOADS=0`, `IRV_METRICS_COUNTRIES=KEN`,
`IRV_METRICS_STATEMENT_TIMEOUT_MS` (default 30000).

## Tests (in pawehe)

```sh
uv run pytest                    # without DATABASE_URL: unit and contract tests; with it (a loaded cube), also the
                                 # /v1 tests, boundary-set and CSV-load tests, and a 1% replay through /v1
                                 # (REPLAY_SAMPLE=1 for all of it)
uv run python -m contract.build --check                                     # openapi.json matches contract/v1.py
uv run python -m tests.replay ../api-data/KEN ../api-data/MOZ --sample 1   # replay the scrape through /v1 (running API)
uv run python -m irv_metrics.check ../api-data/KEN ../api-data/MOZ         # replay against the legacy routes
```
