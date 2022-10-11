# API backend

Server app, written in Python, includes database definition, etl and Tileserver.

Features API leverages PostGres

Tiles API leverages MySQL (until terracotta releases PostGres support officially)

## Installation / Build

`Dockerfile` shows example setup. See
[docs](https://pipenv.pypa.io/en/latest/basics/#pipenv-and-docker-containers)
for how to build venv and run in a more production-oriented way.

Outline of dependencies:
- system `python` and `pip` (e.g. `apt install python3 python3-pip` on Ubuntu)
- `pipenv` or `pip` to manage dependencies (e.g. `pip install pipenv`)
- optionally [`pyenv`](https://github.com/pyenv/pyenv) to provide alternative
  version to system python
- `pyproject.toml` defines python package dependencies
  - includes editable install of `backend` for this application as a package, as
    minimally configured in `setup.py`
- system library dependencies of python packages include
  `libgdal-dev libgeos-dev libpq-dev libproj-dev`


## Configuration

Environment variables:
- use `.env` to define environment variables
- use [`PG*`](https://www.postgresql.org/docs/current/libpq-envars.html) to
  define database connection details. See `.env.example` for an example

```
PYTHONPATH=/code/backend
# Features API
PGDATABASE=jamaica
PGUSER=docker
PGPASSWORD=docker
PGHOST=localhost
# Tiles API
RASTER_BASE_PATH=/data
TC_DRIVER_PATH=mysql://foo:bar@tiles-db/db
TC_SQL_USER=foo
TC_SQL_PASSWORD=bar
TC_ALLOWED_ORIGINS_METADATA='["*"]'
TC_ALLOWED_ORIGINS_TILES='["*"]'
TC_PNG_COMPRESS_LEVEL=0
TC_RESAMPLING_METHOD="nearest"
TC_REPROJECTION_METHOD="nearest"
```

### Tileserver

Tileserver endpoints `/tiles/*` require Cloud-Optimised Tiffs to be mounted at the same path as they have been ingested into the MySQL database.

(Tileserver endpoints require Tiffs to have been pre-loaded into the given MySQL DB)

The base-path for rasters must be set using the environment variable `RASTER_BASE_PATH`

__NOTE__: TC_DRIVER_PATH is not used internally - for Terracotta the path is built programatically based on the URL

Tileserver also provides a meta store for information about each tile database, with associated CRUD operations for metadata management.

#### Adding a Source to the Tileserver metastore

POST the following JSON (adapt as required for the source) to: `http://backend-host:8080/tiles/sources`

##### Aqueduct:

```bash
[
	{
		"source_db": "aqueduct", # the MySQL database the source was ingested-into
		"global_type": "Hazard", # The global hazard type listed for the source tiles
		"domain": "fluvial", # The domain within the UI that the hazard maps-into
		"full_name": "Hazard Aqueduct - Fluvial", # Currently for internal description only
		"description": "description", # Currently for internal description only
		"license": "license", # Currently for internal description only
		"variables": {
			"some": "vars"
		} # Currently for internal description only
	},
	{
		"source_db": "aqueduct",
		"global_type": "Hazard",
		"domain": "coastal",
		"full_name": "Hazard Aqueduct - Coastal",
		"description": "description",
		"license": "license",
		"variables": {
			"some": "vars"
		}
	}
]
```

##### ISIMP Extreme Heat:

```bash
{
	"source_db": "extreme_heat",
	"global_type": "Hazard", 
	"domain": "extreme_heat", 
	"full_name": "Hazard Extreme Heat", 
	"description": "description", 
	"license": "license", 
	"variables": {"some": "vars"} 
}
```