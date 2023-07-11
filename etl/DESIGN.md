# Extract Transform Load (ETL) process design

The contents of this directory provide functionality to process raw data from
various sources into the databases and tile servers that comprise the backend
of infra-risk-vis.

## Datasets

Our datasets can be split into raster and vector types. The majority are
currently raster.

### Raster

Raster metadata is loaded into both the `postgreSQL` database and a `MySQL`
database, tiles are served by the `terracotta` tileserver.

- aqueduct
- esa_land_cover
- exposure_nature
- gem_earthquake
- gem_exposure
- ghsl_buildings
- iris
- isimp_drought
- isimp_extreme_heat
- jrc_pop
- storm
- traveltime_to_healthcare

#### Current approach

For each and every raster dataset:

##### Process raw data into rasters -- snakemake container
- Write an `etl` pipeline to process raw data into cloud optimised geotiffs (COG)
- Adjust `snakemake` service command in docker compose file to target dataset
- Run `snakemake` container to generate COG rasters

##### Raster ingestion -- raster-tile-ingester container
- Bring up the tiles database container, `tiles-db`
- Adjust `raster-tile-ingester` service command in docker compose file to target dataset
- Run `raster-tile-ingester` to use `terracotta` to ingest rasters and create
  metadata in `tiles-db`

##### Metadata creation -- db and backend containers
- Bring up the postgres database container, `db`
- Bring up the backend API container, `backend`
- HTTP POST some high-level raster metadata to the `backend`, in turn creating
  a metadata record in `db`

#### Proposal A

Create a new container for processing, ingestion and metadata creation. The
Dockerfiles for `snakemake` and `raster-tile-ingester` are currently very
similar, with `snakemake` being a geospatial Python environment (with
`terracotta`) and `raster-tile-ingester` being a Python environment with almost
only `terracotta` installed.

Each dataset would have its own service (unique entrypoint). These services
could be collated in a `docker-compose-ingest-rasters.yaml` compose file,
extending the core backend stack from `docker-compose-dev.yaml`, for example.

The entrypoint for each dataset ingester service would be a bash script. The
script will be specific to each dataset as the required processing does vary.
The bash scripts for every dataset would be present in the container.

Each entrypoint bash script will have the following rough outline:
- Calling `snakemake` to run the COG generation (idempotent). Ideally these
  rules would be added to download source data if possible.
- Call `ingest.py` from current `raster-tile-ingester` container to create
  `tiles-db` records.
- Make a POST to `backend` from the new container using the bundled JSON.

#### Proposal B

As in `A`, have a single Dockerfile to produce a container with the
responsibility for the three main stages of raster ingestion. However, extend
`A` by including a layer in the `snakemake` rule paths for dataset, allowing us
to reuse rules across datasets (e.g. cloud optimisation). For those rules which
create records but not files, produce a dummy output to indicate when a DB was
last touched.

The default, `all` rule in this `snakemake` workflow would aim to create files
and DB records for all datasets. The container could take a dataset as an
optional argument (via `docker run`?) to process a single, or set of datasets.

### Vector

Features are loaded into a `postgreSQL` database and served by a vector tileserver.

- osm_rail
- osm_roads
