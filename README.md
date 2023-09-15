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

# Usage

## Data preparation

The visualisation tool runs using prepared versions of analysis data and
results

- Rasters stored as Cloud-Optimised GeoTIFFs, with metadata ingested into
  a terracotta MySQL database, hosted within the backend API.

- Vector data stored in a PostgreSQL database, and preprocessed into Mapbox
  Vector Tiles

See [ETL](etl/README.md) directory for details.

Data to be served from the vector and raster tileservers should be placed on
the host within `tileserver/<data_type>`. These folders are made available to the
running tileservers as docker bind mounts.

For example, in `tileserver/raster/data/<flooding>` there might live TIF files like these:

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

### docker-compose-dev.yaml

Containers used for local development of the entire stack, including ETL.

### docker-compose-prod.yaml

Used to run local builds of Production containers.

### docker-compose-deploy.prod

Used for deployment of containers into a production environment.

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

## Example adding data layer

To add a raster data layer (for example, the `iris` set of tropical cyclone
return period maps):

1. Write an `etl` pipeline, run using `snakemake` container
1. Bring up the tiles database container, `tiles-db`
1. Run `raster-tile-ingester`
1. Bring up the postgres database container, `db`
1. Bring up the backend API container, `backend`
1. Post the raster metadata to the backend - see [docs](containers/backend/README.md)
   on adding data to the tileserver metadata store.

## Example service update

Locally, build the frontend, push to the container registry:

```bash
# Edit docker-compose.yaml image version, in this example line 28:
#     image: ghcr.io/nismod/gri-web-server:0.16

# Build
docker compose -f docker-compose-prod.yaml build web-server

# Log in to the container registry
# see: https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry

# Push
docker push ghcr.io/nismod/gri-web-server:0.16
```

On remote, pull the image and reload service:

```bash
# Pull image
docker pull ghcr.io/nismod/gri-web-server:0.16

# Edit docker-compose.yaml image version (or sync up), in this example line 28:
#     image: ghcr.io/nismod/gri-web-server:0.16

# Restart service
docker compose up -d web-server
```

## IRV AutoPackage Service

Provides API for extraction of data (and hosting of results) from various layers using pre-defined boundaries.

See [`irv-autopkg`](http://github.com/nismod/irv-autopkg) for more information.

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
  and Global Earthquake Model Foundation. This work has been funded by the World
  Bank Group, Willis Towers Watson Insurance for Development Forum, the UK
  Natural Environment Research Council (NERC) through the UK Centre for Greening
  Finance and Investment, and UK Foreign, Commonwealth and Development Office
  (FCDO) through the Climate Compatible Growth (CCG) programme.
