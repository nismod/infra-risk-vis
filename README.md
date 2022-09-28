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
  a terracotta SQLite database
- Vector data stored in a PostgreSQL database, and preprocessed into Mapbox
  Vector Tiles

See `./etl` directory for details.

Data to be served from the vector and raster tileservers should be placed on
the host within `mounts/<data_type>`. These folders are made available to the
running tileservers as docker bind mounts.

For example, in `mounts/raster/` there might live TIF files like these:
```
coastal_mangrove__rp_100__rcp_baseline__epoch_2010__conf_None.tif
coastal_mangrove__rp_25__rcp_baseline__epoch_2010__conf_None.tif
coastal_mangrove__rp_500__rcp_baseline__epoch_2010__conf_None.tif
coastal_nomangrove_minus_mangrove__rp_100__rcp_baseline__epoch_2010__conf_None.tif
```

And in `mounts/vector/`, mbtiles files like these:
```
airport_runways.mbtiles
airport_terminals.mbtiles
buildings_commercial.mbtiles
buildings_industrial.mbtiles
```

### Data preperation within Docker

Data Preperation can be run within Docker:

bash
```
docker run -it -v ${PWD}/etl:/opt/etl gri-snakemake:latest --cores 1 -s /opt/etl/Snakefile
```

or using Docker Compose `snakemake` service:

bash
```
docker-compose -f docker-compose-dev.yaml run snakemake
```

### Terracotta Test Frontend

Boot the `raster-tileserver` service then run:

```bash
terracotta connect localhost:5000
```

## Build

The application is built with several 'services', each facilitated by a running
docker container.

Services:
- Web server (nginx)
- Raster tileserver (terracotta)
- Vector tileserver (tileserver-gl)
- Backend / API (bespoke Python app)
- Database (PostgreSQL with PostGIS)

The services are orchestrated using `docker compose`. N.B. The app was built
with docker engine version 20.10.16 and compose version 2.5.0. It may not work
with other versions.

The `compose.yml` file contains service names and definitions, and paths to
build contexts (all located within `containers/`).

To build images for all the services, use `docker compose build`, and to build
for a single service, use `docker compose build <service>`, e.g.
`docker compose build raster-tileserver`. These commands are automatically
referring to the `compose.yml` file, meaning `docker compose build` is
equivalent to `docker compose -f compose.yml build`.

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

The application has access control. Authentication is handled by nginx. To add
a user, first exec into the web-server container. The following command will
start a `/bin/sh` shell in that container:
`docker exec -it $(docker ps | grep web-server | awk '{printf $1}') /bin/sh`

Then create a new user and password pair with `htpasswd`:
`htpasswd -Bb /etc/nginx/auth/.htpasswd <username> <password>`

The passwords are persisted as a docker volume. If this volume is deleted, the
accounts will be lost.

The site can run on a single Linux machine or virtual machine. See `./deploy`
directory for further details.

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
