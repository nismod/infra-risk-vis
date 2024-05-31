# Extract Transform Load (ETL) for infra-risk-vis (IRV)

## Introduction

This directory contains the logic, configuration and documentation for
acquiring, processing and ingesting data such that a running IRV stack may host
as tilesets.

The data processing steps are broadly as follows:

### Raster

- Download data if possible
- Reproject to WGS84 if necessary
- Set zeros to no data value
- Clip to remove polar regions
- Cloud optimise
- Ingest into `terracotta`, creating a metadata database for dataset if
  necessary
- Create dataset metadata in Postgres database

### Vector

Vector data are yet to be incorporated in the unified ETL workflow. See the
relevant pipeline readme files for more information.

## Architecture

The ETL process is currently implemented with a combination of scripts and
[snakemake](https://snakemake.readthedocs.io/en/stable/) rules. There is
currently an effort underway to move _all_ processing logic into snakemake
rules.

The `Snakefile` and `.smk` files contain the rule definitions for deciding how
to transform some input file into an output. Many datasets have some of their
own specific rules, typically for downloading and initial processing, for
instance, reprojection.

## Setup

### Software environment

Install the necessary dependencies:

```bash
micromamba create -f environment.yml -y
```

To activate the environment:

```bash
micromamba activate irv-etl
```

### Required services

The last two stages of the ETL pipeline (ingestion and metadata creation)
involve interacting with the database service defined in the parent directory to
this one.

To bring up these services, refer to the [readme](../README.md) in the parent
directory for a full explanation of their required env files, etc., but briefly:

```bash
docker compose -f docker-compose-dev.yaml up db -d
```

### Awkward files

Unfortunately, not all the source data is openly available on the internet. Some
files must be manually copied into the appropriate location prior to running the
ETL process.

The affected datasets include:

- gem_earthquake
- iris

Check the `snakemake` rules for more information, but source raster data should
typically reside in `raster/raw/<dataset>/`.

You may wish to remove write permissions to these files once they have been
installed, e.g. `chmod ug-w raster/raw/gem_earthquake/*.tif`. This means `rm -r
raster` will remove files than can be replaced automatically, but not the
awkward files.

### Configuration

The ETL pipeline is primarily configured with an environment file, located at
`../envs/{dev|prod}/.etl.env`. This contains details for connections to services
and authentication information for pipelines which require it.

Here is an example environment file to use as a template:

```bash
# configuration for terracotta
TC_DRIVER_PATH=postgresql://global_dev:password@localhost:5432
TC_DRIVER_PROVIDER=postgresql
TC_PNG_COMPRESS_LEVEL=0
TC_RESAMPLING_METHOD=nearest
TC_REPROJECTION_METHOD=nearest

# connecting to database
PGHOST=localhost
PGDATABASE=global_dev
PGUSER=global_dev
PGPASSWORD=password

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
snakemake --cores <n_cores> -- raster/cog/<dataset>/<key>.tif
```

For example:

```bash
snakemake --cores <n_cores> -- raster/cog/exposure_nature/ocs_0-30cm_mean_1000.tif
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
snakemake --cores <n_cores> -- raster/metadata/<dataset_name>.flag
```

### Running every pipeline

To run every pipeline, we do not request a file, but rather a target rule called
`all`.

```bash
snakemake --cores <n_cores> -R all
```

This will create and ingest all the pertinent rasters and create metadata
records for them.

## Adding new datasets

### Raster

To add additional raster datasets to the ETL pipeline you will need several new
files in `pipelines/<new_dataset_name>`:

- `README.md` describing the dataset
- `layers.csv` table of layers
- `metadata.json` metadata to be written to the Postgres database
- `rules.smk` to acquire and process the data into a WGS84 raster. The rule(s)
  should write output files to `raster/raw/<dataset>/`. The shared rules will
  then clip, set zero to no data, and cloud optimise unless overridden by custom
  rules.

The `layers.csv` file must follow the example structure:

- `filename` column must contain the file basename (no directory path)
- one or more columns, to be listed as `keys` in the `metadata.json` must exist
  and be sufficient to uniquely identify each row
- additional columns may exist (e.g. URL) to support the raster workflow but
  will be ignored in the metadata processing

```csv
filename,mode
202001_Global_Motorized_Travel_Time_to_Healthcare_2019.tif,motorized
202001_Global_Walking_Only_Travel_Time_To_Healthcare_2019.tif,walking
```

The `metadata.json` must follow the example structure:

- `"domain"` must match the `DATASET` wildcard used in the the snakemake rules
- `"name"` should be a short, readable description
- `"group"` should be one of the high-level front-end groups (Hazard, Exposure,
  Vulnerability, Risk)
- `"description"` should be a short description or citation
- `"license"` should be a short code for the open data license
- `"keys"` must be the list of metadata keys to be loaded to the terracotta
  metadata database, identical to the variable column names used in `layers.csv`
  to identify/address each layer

```json
{
  "domain": "traveltime_to_healthcare",
  "name": "Global maps of travel time to healthcare facilities",
  "group": "Exposure",
  "description": "Weiss, D.J., Nelson, A., Vargas-Ruiz, C.A. et al. Global maps of travel time to healthcare facilities. Nat Med 26, 1835–1838 (2020). https://doi.org/10.1038/s41591-020-1059-1",
  "license": "CC-BY 4.0",
  "keys": ["mode"]
}
```

The `Snakefile` will also require modification:

- If you have written a new `rules.smk`, you will need to import it as a module
  here. See existing datasets for more information.
- If you wish to overwrite the behaviour of a common rule, e.g. clipping, cloud
  optimsation, etc, you can override rules when importing using `snakemake`'s
  `ruleorder` directive.
- You should also add your dataset to `ALL_DATASETS` so that the `all` target
  rule will work as expected.

### Vector

To be implemented and then documented!
