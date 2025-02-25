# Infrastructure Risk Visualisation Tool

This project provides interactive data visualisations of risk analysis results.

![About](images/screenshot-about.png)

The tool presents the infrastructure systems and hazards considered in the
analysis, then presents results as modelled for the whole system at a fine
scale.

Other functionality:

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

# Architecture

The tool runs as a set of containerised services:

- Traefik reverse proxy to direct requests to the other services
- Web server (nginx) `ghcr.io/nismod/gri-web-server`
- Vector tileserver (tileserver-gl) `ghcr.io/nismod/gri-vector-tileserver`
- Backend / API (bespoke Python app for vector data and raster tiles (+meta)) `ghcr.io/nismod/gri-backend`
- API Database (Postgres with PostGIS serves backend) (Dev only)
- Tiles Database (Postgres server with multiples terracotta metadata databases) (Dev only)

The services are orchestrated using `docker compose`.

N.B. The app was built with docker engine version 20.10.16 and compose version
2.5.0. It may not work with other versions.

# Usage

## Data preparation

The visualisation tool runs using prepared versions of analysis data and
results:

- Rasters stored as Cloud-Optimised GeoTIFFs, with metadata ingested into
  a terracotta database, hosted within the backend API.
- Vector data stored in a PostgreSQL database, and preprocessed into Mapbox
  Vector Tiles

See [ETL](etl/README.md) directory for details.

Data to be served from the vector and raster tileservers should be placed on
the host within `tileserver/<data_type>`. These folders are made available to the
running tileservers as docker bind mounts.

For example, in `tileserver/raster/data/aqueduct` there might live TIF files like these:

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

## Environment

Environment variables for the various services (and the ETL workflow) are
stored in env files. Example files are given in `envs/dev-example`. These can
be placed in `envs/dev` to get started.

Production env files should be placed in `envs/prod`.

## Deploy

To deploy the stack we use the `docker compose` tool.

### Development

The set of long-running services can include:

- traefik: Reverse proxy for other services, handles TLS
- web-server: Nginx server for the frontend code and static files
- db: PostgreSQL database holding vector data and raster metadata
- backend: API for available datasets and raster tileserver (terracotta)
- vector-tileserver: TileServer-GL for serving .mbtiles files
- redis: In-memory database for autopackage job queueing
- irv-autopkg-worker: Autopackage data processing (clipping, serialisation, etc.)
- irv-autopkg-api: Autopackage service coordination

If you're running your own [frontend](https://github.com/nismod/irv-frontend/)
development server, or connecting to a remotely hosted database, or not using
the [autopackage API](https://github.com/nismod/irv-autopkg), you may not need
all these services.

To this end, we use [profiles](https://docs.docker.com/compose/profiles/) to
define 'core' services which always run, and optional services. A bare `docker
compose -f docker-compose-dev.yaml up` will run only the core services (those
without a `profiles` attribute).

For example, when running your own FE development server to add a new raster
layer the following should suffice: `docker compose -f docker-compose-dev.yaml
up`. This will bring up `db`, `tiles-db`, `backend` and `vector-tileserver`.

To run the core services with a standard frontend: `docker compose -f
docker-compose-dev.yaml --profile web-server up`.

To run the core services alongside the autopackage services: `docker compose -f
docker-compose-dev.yaml --profile autopkg up`.

To run all of these behind traefik (every long-running service): `docker
compose -f docker-compose-dev.yaml --profile traefik --profile web-server
--profile autopkg up`.

There are also a few short-lived 'utility containers', which can be run to
perform particular tasks:

- recreate-metadata-schema: Drop the contents of the `db` database, recreate with empty tables
- raster-tile-delete-entries: Delete raster entries of specified dataset in `tiles-db`
- raster-tile-drop-database: Drop whole database for specified dataset from `tiles-db`

When starting from a clean slate, the `recreate-metadata-schema` service must
be run to create the tables in `db` that `backend` relies upon. If you find
that the `backend` service is complaining that the `raster_tile_sources`
database table is not available, you may need to create the appropriate tables
in the `db` service first. To do that, bring the `db` service up as described
above, and then run: `docker-compose -f docker-compose-dev.yaml up
recreate-metadata-schema` to (re)create the tables. Note that this will drop
any data currently in database.

### Production

To run local builds of production containers we use the
`docker-compose-prod-build.yaml` file. See [below](#Updating a service) for
more details.

To deploy containers into a production environment:
`docker compose -f docker-compose-prod-deploy.yaml up -d`

## Updating a service

To update a service:

- We make the necessary changes to the container
- Build a new container
- Push it to the container repository
- Pull it on the production machine
- Deploy it

As an example, below we update the backend on a development machine:

```bash
# Edit docker-compose-prod-build.yaml image version:
#     image: ghcr.io/nismod/gri-backend:1.8.1

# Build
docker compose -f docker-compose-prod-build.yaml build backend

# Log in to the container registry
# see: https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry

# Push
docker push ghcr.io/nismod/gri-backend:1.8.1
```

On the production remote, pull the image and restart the service:

```bash
# Pull image
docker pull ghcr.io/nismod/gri-backend:1.8.1

# Edit docker-compose-prod-deploy.yaml image version (or sync up):
#     image: ghcr.io/nismod/gri-backend:1.8.1

# Restart service
docker compose -f docker-compose-prod-deploy.yaml up -d backend
```

## Adding new data layers

To add a raster data layer (for example, the `iris` set of tropical cyclone
return period maps) see the [ETL](etl/README.md) directory.

## IRV AutoPackage Service

Provides API for extraction of data (and hosting of results) from various
layers using pre-defined boundaries.

See [`irv-autopkg`](http://github.com/nismod/irv-autopkg) for more information.

# Acknowledgements

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
- [v0.3](https://github.com/oi-analytics/oi-risk-vis/releases/tag/v0.3.0-jamaica)
  was developed by the Oxford Programme for Sustainable Infrastructure Systems
  (OPSIS) in the Environmental Change Institute, University of Oxford, for the
  Government of Jamaica (GoJ) as part of a project funded by UK Aid (FCDO). The
  initiative forms part of the Coalition for Climate Resilient Investment’s
  (CCRI) collaboration with the GoJ, which also includes analysis of
  nature-based approaches to build resilience in Jamaica to be procured and
  funded by the Green Climate Fund (GCF).
- [release/caribbean](https://github.com/nismod/infra-risk-vis/tree/release/caribbean)
  was developed as part of the Jamaica project.
- [release/east-africa](https://github.com/nismod/infra-risk-vis/tree/release/east-africa)
  was developed by researchers in the University of Southampton's Transportation
  Research Group and the Oxford Programme for Sustainable Infrastructure
  Systems, University of Oxford, supported by engagement with infrastructure and
  climate specialists and related government bodies, and funded by UKAID through
  the UK Foreign, Commonwealth & Development Office under the High Volume
  Transport Applied Research Programme, managed by DT Global.
- current work on global-scale data and analysis continues to be led by
  researchers in OPSIS. In part this is in collaboration with the Global
  Resilience Index Initiative, including the Oxford Spatial Finance Initiative
  and Global Earthquake Model Foundation. This work has been funded by: the World
  Bank Group; Willis Towers Watson Insurance for Development Forum; the UK
  Natural Environment Research Council (NERC) through the UK Centre for Greening
  Finance and Investment; the UK Foreign, Commonwealth and Development Office
  (FCDO) through the Climate Compatible Growth (CCG) programme; Howden Foundation;
  and Global Center on Adaptation.
