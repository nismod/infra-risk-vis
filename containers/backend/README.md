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
cd containers/backend

# create a virtual environment (using venv or another method if you prefer)
python -m venv venv
source venv/activate
pip install -e .[dev]

# run the application
uvicorn app.main:app --port 8888 --reload
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
