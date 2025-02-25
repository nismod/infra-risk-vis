# GRI Metrics API

REST API, written in Python, includes database definition.

The API is a light wrapper over tables stored in Postgres, using the
PostGIS extension for spatial data types.

## Installation / Build

### Docker

```bash
# Build the metrics service
docker-compose -f docker-compose-prod-build.yaml build metrics-api
# Start the database to run in the background
docker-compose -f docker-compose-dev.yaml up -d database
# Start the metrics api to run in the foreground - CTRL+C to quit
docker-compose -f docker-compose-dev.yaml up metrics-api
```

### Running Locally (alternative to Docker)

```bash
# cd to this backend app directory
cd containers/metrics

# create a virtual environment (using venv or another method if you prefer)
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# run the application (from /metrics)
uvicorn api.main:app --port 8888 --reload
```

## Configuration

Environment variables:

- use `.env` to define environment variables
- use [`PG*`](https://www.postgresql.org/docs/current/libpq-envars.html) to
  define database connection details. See `.metrics.env` for an example.

```
METRICS_DB_URL=
METRICS_LOG_LEVEL=
```

## Datasets

### GDL

Data downloaded from Global Data Lab

API base path: /metrics/gdl

subnational GeoJson downloaded/simplified from (V6.4) shapefiles:
https://globaldatalab.org/mygdl/downloads/shapefiles/

national GeoJson downloaded from visualisation:
https://globaldatalab.org/geos/maps/

annual development, education, income, and healthcare metrics downloaded from:
https://globaldatalab.org/shdi/download/

ETL scripts and data are provided in etl/gdl.
