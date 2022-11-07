# API backend

Server app, written in Python, includes database definition, etl and Tileserver.

Features API leverages PostGres

Tiles API leverages Terracotta Python API and MySQL.

Tiles must be loaded sperately - there are no endpoints for ingesting data at present.

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
PGDATABASE=
PGUSER=
PGPASSWORD=
PGHOST=

# Tiles API
RASTER_BASE_PATH=/data # Path at-which raster tiles can be found (must match the MySQL Tiles-db loaded path)
TC_DRIVER_PATH=mysql://foo:bar@tiles-db/db # Tiles-db MySQL Host
TC_SQL_USER=foo
TC_SQL_PASSWORD=bar
TC_ALLOWED_ORIGINS_METADATA='["*"]'
TC_ALLOWED_ORIGINS_TILES='["*"]'
TC_PNG_COMPRESS_LEVEL=0
TC_RESAMPLING_METHOD="nearest"
TC_REPROJECTION_METHOD="nearest"

API_TOKEN= # API token is only required for mutation operation.
DOMAIN_TO_DB_MAP='{\"land_cover\":\"land_cover\"}' # Valid JSON of a mapping between front-end DOMAIN values and the database in-which the data is stored.
```

### Tileserver

Tileserver endpoints `/tiles/*` require Cloud-Optimised Tiffs to be mounted at the same path as they have been ingested into the MySQL database.

(Tileserver endpoints require Tiffs to have been pre-loaded into the given MySQL DB)

The base-path for rasters must be set using the environment variable `RASTER_BASE_PATH`

__NOTE__: TC_DRIVER_PATH is not used internally - for Terracotta the path is built programatically based on the URL

Tileserver also provides a meta store for information about each tile database, with associated CRUD operations for metadata management.

#### Categorical Data

Categorical rasters are supported.  Categorical colormaps can either be included in the `config.py`, or passed with each tile request, as per the terracotta documentation:  https://terracotta-python.readthedocs.io/en/latest/tutorials/categorical.html

__NOTE__: Only `{pixel :(RGBA)}` explicit color maps are supported in either the request or `config,py`

If included in `config.py` the key must match an existing MySQL database, with pre-loaded categorical raster(s).

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

__NOTE__: Large rasters can fail to load due to dropped MySQL connections to Cloud-hosts after metadata creation.  This appears to be a bug with the underlying library.  The fix to-date has been to load them locally and subsequently push the MySQL database to the Cloud host.

#### Adding a Source to the Tileserver metastore

1. Add the source to the domain:mysqldatabase mapping in the config (which is loaded from the environment as json):

```json
{
    "fluvial": "aqueduct",
    "coastal": "aqueduct",
    "extreme_heat": "extreme_heat",
    "cyclone": "cyclone",
}
```

E.g. in the environment one would add the following:

```
DOMAIN_TO_DB_MAP='{\"land_cover\":\"land_cover\",...}'
```

2. `POST` metadata about the raster to: `http://backend-host:8080/tiles/sources`.  Examples of the payload are show below.


```json
{
  "source_db": "aqueduct", # the MySQL database the source was ingested-into
  "global_type": "Hazard", # The global hazard type listed for the source tiles
  "domain": "fluvial", # The domain within the UI that the hazard maps-into
  "full_name": "Hazard Aqueduct - Fluvial", # Currently for internal description only
  "description": "description", # Currently for internal description only
  "license": "license", # Currently for internal description only
  "variables": {} # Currently for internal description only
}
```

##### Aqueduct:

```json
[
	{
		"source_db": "aqueduct",
		"global_type": "Hazard",
		"domain": "fluvial",
		"full_name": "Hazard Aqueduct - Fluvial",
		"description": "description",
		"license": "license",
		"variables": {}
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

```json
{
	"source_db": "extreme_heat",
	"global_type": "Hazard", 
	"domain": "extreme_heat", 
	"full_name": "Hazard Extreme Heat", 
	"description": "description", 
	"license": "license", 
	"variables": {
		"gcm": "gcm",
		"rcp": "rcp",
		"type": "hazard",
		"epoch": "epoch",
		"metric": "metric"
	}
}
```

##### STORM

```json
{
	"source_db": "storm",
	"global_type": "Hazard", 
	"domain": "cyclone", 
	"full_name": "Hazard Tropical Storm", 
	"description": "description", 
	"license": "license", 
	"variables": {
		"type": "hazard",
		"rp": "rp",
		"gcm": "gcm"
	}
}
```

##### JRC Population

```json
{
  "source_db": "jrc_pop",
  "global_type": "Exposure",
  "domain": "population",
  "full_name": "JRC Population",
  "description": "description",
  "license": "license",
  "variables": {}
}
```

##### ISIMP Drought

```json
{
  "source_db": "drought",
  "global_type": "Hazard",
  "domain": "drought",
  "full_name": "ISIMP Drought",
  "description": "description",
  "license": "license",
  "variables": {}
}
```

##### Exposure Nature

```json
{
  "source_db": "exposure_nature",
  "global_type": "Exposure",
  "domain": "nature",
  "full_name": "Nature Exposure",
  "description": "description",
  "license": "license",
  "variables": {}
}
```

##### GDSL Buildings

```json
{
  "source_db": "buildings",
  "global_type": "Exposure",
  "domain": "buildings",
  "full_name": "Building Exposure",
  "description": "description",
  "license": "license",
  "variables": {}
}
```

##### Traveltime To Healthcare

```json
{
  "source_db": "traveltime_to_healthcare",
  "global_type": "Exposure",
  "domain": "traveltime_to_healthcare",
  "full_name": "Travel Time to Healthcare",
  "description": "description",
  "license": "license",
  "variables": {}
}
```

##### ESA Land Cover

```json
{
  "source_db": "land_cover",
  "global_type": "Exposure",
  "domain": "land_cover",
  "full_name": "ESA Land Cover",
  "description": "description",
  "license": "license",
  "variables": {}
}
```