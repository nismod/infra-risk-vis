# Extract Transform Load (ETL) for infra-risk-vis (IRV)

## Introduction

This directory contains the logic, configuration and hopefully documentation for
acquiring, processing and ingesting data such that a running IRV stack may host
as tilesets.

The data processing steps are broadly as follows:

### Raster

- Download data if possible
- Reproject to WGS84 if necessary
- Set zeros to no data value
- Clip to remove polar regions
- Cloud optimise
- Ingest into `terracotta`, creating `mysql` database for dataset
- Create dataset metadata in `postgreSQL` database

### Vector

Vector data are yet to be incorporated in the unified ETL workflow. See the
relevant pipeline readme files for more information.

## Architecture

The ETL process is currently implemented with a combination of scripts and
[snakemake](https://snakemake.readthedocs.io/en/stable/) rules. There is currently an
effort underway to move _all_ processing logic into snakemake rules.

The `Snakefile` and `.smk` files contain the rule definitions for deciding how
to transform some input file into an output. Many datasets have some of their
own specific rules, typically for downloading.

## Setup

### Software environment

Install the necessary dependencies using conda/mamba/micromamba.

```bash
micromamba create -f environment.yml -y
```

To activate the environment:

```bash
micromamba activate irv-etl
```

### Required services

The later stages of the ETL pipeline involve interacting with services
defined in the parent directory to this one. The `backend` service acts as an
intermediary to two databases, a `postgreSQL` and a `mysql` instance, known as
the `db` and `tiles-db` services respectively.

To bring up these services, refer to the [readme](../README.md) in the parent
directory for a full explanation of their required env files, etc., but briefly:

```bash
docker compose -f docker-compose-dev.yaml up db tiles-db backend
```

### Awkward files

Unfortunately, not all the source data is openly available on the internet.
Some files must be manually copied into the appropriate location prior to
running the ETL process.

The affected datasets include:
- gem_earthquake
- iris

Check the `snakemake` rules for more information, but source raster data should
typically reside in `raster/raw/<dataset>/`.

You may wish to remove write permissions to these files once they have been installed.

### Configuration

The ETL pipeline is primarily configured with an environment file, located at
`../envs/{dev|prod}/.etl.env`. This contains details for connections to services
and authentication information for pipelines which require it.

Here is an example environment file to use as a template:

```bash
# connecting to `tiles-db`, MySQL
TC_DRIVER_PATH=mysql://root:password@tiles-db
TC_DRIVER_PROVIDER=mysql
TC_PNG_COMPRESS_LEVEL=0
TC_RESAMPLING_METHOD=nearest
TC_REPROJECTION_METHOD=nearest

# connecting to `db`, postgreSQL
PGHOST=db
PGDATABASE=global
PGUSER=docker
PGPASSWORD=docker

# connecting to `backend`
BE_HOST=localhost
BE_PORT=8888
BE_API_TOKEN=test  # required for mutation operations on tile metadata (`/tiles/sources POST & DELETE`).

# data downloading
# https://cds.climate.copernicus.eu/api-how-to
COPERNICUS_CDS_URL="https://cds.climate.copernicus.eu/api/v2"
COPERNICUS_CDS_API_KEY=  # "<uid>:<token>"
```

## Usage

With the software environment activated (see above), one can request files via
`snakemake` to invoke jobs. These may be processed rasters, or in the case of
database operations, dummy files with the extension `.flag`. Requesting any
missing output implies its ancestors are also required.

### Running a single raster

To request a cloud optimised raster, invoke `snakemake` as follows:

```bash
snakemake --cores <n_cores> -- rasters/cog/<dataset>/<key>.tif
```

For example:

```bash
snakemake --cores <n_cores> -- rasters/cog/exposure_nature/ocs_0-30cm_mean_1000.tif
```

N.B. The `--dry-run` or `-n` option can be used to preview which jobs
`snakemake` has determined are necessary prior to executing.

To request a raster processed to another stage in the pipeline, substitute
`raw`, `no_data` or `clip` for `cog` in the above paths.

### Running all rasters in a single pipeline

A list of datasets currently implemented in the unified workflow is kept as
`ALL_DATASETS` in the `Snakefile`.

The full processing pipeline for a single raster dataset will acquire and
process the rasters, ingest them and create a metadata record. Invoke it as
follows:

```bash
snakemake --cores <n_cores> -- pipelines/<dataset_name>/metadata_created.flag
```

### Running every pipeline

To run every pipeline, we do not request a file, but rather a target rule called `all`.

```bash
snakemake --cores <n_cores> -R all
```

This will create and ingest all the pertinent rasters and create metadata records for them.

## Extension - adding a new dataset

### Raster

To add additional raster datasets to the ETL pipeline you will need several new
files in `pipelines/<new_dataset_name>`:

- A `README.md` containing describing what the dataset is and where it was
sourced from.
- A `layers.csv` file containing CSV table of layers.  This should include a
header with the `path` of the layer, the `type` and any variables you wish to be
able to filter by.
- A `metadata.json` file containing metadata to be written to the `postgreSQL`
database. The mapping referenced by the `variable` key will be used to configure
`terracotta`'s `mysql` database, including which keys are available to filter
rasters by.
- A `rules.smk` rules file containing rules to acquire and process the data into
a WGS84 raster. The rule(s) should write output files to
`raster/raw/<new_dataset>/`.

The `Snakefile` will also require modification:
- If you have written new rules, you will need to import them as a module here.
See existing datasets for more information.
- If you wish to overwrite the behaviour of a common rule, e.g. clipping, cloud
optimsation, etc, you can override rules when importing.
- You should also add your dataset to `ALL_DATASETS` so that the `all` target
rule will work as expected.

### Vector

To be implemented and then documented!