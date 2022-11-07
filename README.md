# Infrastructure Risk Visualisation Tool

This project provides interactive data visualisations of risk analysis results.

![About](images/screenshot-about.png)

The tool presents the infrastructure systems and hazards considered in the
analysis, then presents results as modelled for the whole system at a fine
scale.

See an overview of infrastructure networks:

![Networks](images/screenshot-overview.png)

Other functionality planned (and incorporated in some way in previous versions):

- Summarise risk analysis at an administrative regional scale.
- Zoom in to see networks in detail.
- See an overview of hazard data.
- Inspect details of hazard layers.
- Query attributes of elements of the system.
- Range of potential economic impacts of failure, consisting of direct damages
  to infrastructure assets and indirect economic losses resulting from
  infrastructure service disruption (loss of power, loss of access).
- Explore a cost-benefit analysis (under uncertainty, with options to explore
  some parameters) of adaptation measures.

This README covers requirements and steps through how to prepare data for
visualisation and how to run the tool.

# Usage

## Data preparation

The visualisation tool runs using prepared versions of analysis data and
results
- Rasters stored as Cloud-Optimised GeoTIFFs, with metadata ingested into
  a terracotta MySQL database, hosted within the backend API.
- Vector data stored in a PostgreSQL database, and preprocessed into Mapbox
  Vector Tiles

See `./etl` directory for details.

Data to be served from the vector and raster tileservers should be placed on
the host within `tileserver/<data_type>`. These folders are made available to the
running tileservers as docker bind mounts.

For example, in `tileserver/raster/` there might live TIF files like these:
```
coastal_mangrove__rp_100__rcp_baseline__epoch_2010__conf_None.tif
coastal_mangrove__rp_25__rcp_baseline__epoch_2010__conf_None.tif
coastal_mangrove__rp_500__rcp_baseline__epoch_2010__conf_None.tif
coastal_nomangrove_minus_mangrove__rp_100__rcp_baseline__epoch_2010__conf_None.tif
```

And in `tileserver/vector/`, mbtiles files like these:
```
airport_runways.mbtiles
airport_terminals.mbtiles
buildings_commercial.mbtiles
buildings_industrial.mbtiles
```

### Docker Development Environment

`docker-compse-dev.yaml` includes a set of services for use in the dataprocessing and development process.

The following environment files are required:

#### Environment

The following env files are required (in `envs/dev/.*`):

##### .backend.env

```
PGHOST=
PGDATABASE=
PGUSER=
PGPASSWORD=

# Tiles API
LOG_LEVEL=INFO
RASTER_BASE_PATH=/data  # The mount underwich GeoTiffs for the tileserver can be found
MYSQL_URI=  # MySQL URI for tiles-db
API_TOKEN=  # Only required for mutating tiles metadata in the API

# Terracotta internal
TC_ALLOWED_ORIGINS_METADATA='["*"]'
TC_ALLOWED_ORIGINS_TILES='["*"]'
TC_PNG_COMPRESS_LEVEL=0
TC_RESAMPLING_METHOD="nearest"
TC_REPROJECTION_METHOD="nearest"
```

##### .db.env

```
POSTGRES_DB=
POSTGRES_USER=
POSTGRES_PASS=
ALLOW_IP_RANGE=0.0.0.0/0
POSTGRES_MULTIPLE_EXTENSIONS=postgis
```

##### .mysql.env

```
MYSQL_USER=
MYSQL_PASSWORD=
MYSQL_ROOT_PASSWORD=
```

##### .pgadmin.env

```
PGADMIN_DEFAULT_EMAIL=
PGADMIN_DEFAULT_PASSWORD=
WORKERS=1
```

##### .raster-tile-ingester.env

```
# Terracotta Env
TC_DRIVER_PATH=mysql://
TC_DRIVER_PROVIDER=mysql
TC_PNG_COMPRESS_LEVEL=0
TC_RESAMPLING_METHOD="nearest"
TC_REPROJECTION_METHOD="nearest"

# Gri Backend Env - for managing entries in the internal API tileserver
BACKEND_HOST=
BACKEND_PORT=
API_TOKEN=
```

### Data preperation within Docker

Data Preperation can be run within Docker, end to end.

pipelines/{workflow}/{workflow_preprocessor.py} (pre-processing and generation of csv for snakemake) -> pipelines/{workflow}/Snakemake (generated COG Files) -> docker raster-tile-ingester (ingests to tile server DB) -> Add meta to backend API (see: containers/backend/README.md)

#### Snakemake

bash
```
docker run -it -v ${PWD}/etl:/opt/etl gri-snakemake:latest --cores 1 -s /opt/etl/Snakefile
```

or using Docker Compose `snakemake` service:

bash
```
docker-compose -f docker-compose-dev.yaml run snakemake
```

#### Raster Tileserver Ingester

Update docker-compose-dev.yaml `raster-tile-ingester` block as req. for a dataset (after running its pre-processing and Snakemake pipeline)

```bash
docker-compose -f docker-compose-dev.yaml run raster-tile-ingester
```


## Build

The application is built with several 'services', each facilitated by a running
docker container.

Services:
- Web server (nginx)
- Vector tileserver (tileserver-gl)
- Backend / API (bespoke Python app for vector data and raster tiles (+meta))
- API Database (PostgreSQL with PostGIS serves backend)
- Tiles Database (MySQL serves tile ingester and backend /tiles endpoints)

The services are orchestrated using `docker compose`. N.B. The app was built
with docker engine version 20.10.16 and compose version 2.5.0. It may not work
with other versions.

The `compose.yml` file contains service names and definitions, and paths to
build contexts (all located within `containers/`).


## Deploy

To deploy the stack, use `docker compose up`. Again, this will default to using
the `compose.yml` file. The current process will show the interleaved log
output of the various services. To bring up the stack in the background (e.g.
for production), use `docker compose up -d` to daemonise the process.

The default service configuration can be modified by means of 'overlaying'
compose files. For development purposes it can be useful to 'bind mount' source
code files into running containers. Similarly, direct access to the containers'
various open ports can help debugging. To deploy the stack locally with such
changes, use the following invocation:
`docker compose -f compose.yml -f compose-local.yml up`

To stop a foregrounded compose stack, issue SIGTERM with Ctrl-C. If services
haven't stopped in 10 seconds they will be brutally terminated. To bring a
daemonised stack down, use `docker compose down`.

## Acknowledgements

This tool has been developed through several projects.

- [v0.1](https://github.com/oi-analytics/oi-risk-vis/releases/tag/v0.1-argentina)
  was developed by Oxford Infrastructure Analytics for the Government of
  Argentina with funding support from the World Bank Group and Global Facility
  for Disaster Reduction and Recovery (GFDRR).
- [v0.2](https://github.com/oi-analytics/oi-risk-vis/releases/tag/v0.2.0-seasia)
  was developed by Oxford Infrastructure Analytics for the Disaster Risk
  Financing and Insurance Program (DRFIP) of the World Bank with support from
  the Japan&mdash;World Bank Program for Mainstreaming DRM in Developing
  Countries, which is financed by the Government of Japan and managed by the
  Global Facility for Disaster Reduction and Recovery (GFDRR) through the Tokyo
  Disaster Risk Management Hub.
- current development is by the Oxford Programme for Sustainable Infrastructure
  Systems in the Environmental Change Institute, University of Oxford, for the
  Government of Jamaica (GoJ) as part of a project funded by UK Aid (FCDO). The
  initiative forms part of the Coalition for Climate Resilient Investment’s
  (CCRI) collaboration with the GoJ, which also includes analysis of
  nature-based approaches to build resilience in Jamaica to be procured and
  funded by the Green Climate Fund (GCF).
