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
relevant pipeline readme files for more information. The rest of this section
describes the general approach.

We use [tippecanoe](https://github.com/mapbox/tippecanoe) to generate Mapbox Vector
Tiles, stored in `.mbtiles` files. Follow the installation or build instructions
in their documentation.

Step 1 is to have the features prepared in GeoJSON.

The simplest `tippecanoe` example command often works well enough:

```bash
tippecanoe -zg -o landslide_forest.mbtiles --drop-densest-as-needed landslide_forest.geojson
```

Here's a version with options that should work a little better for a larger dataset:

```bash
tippecanoe  -o landslide_forest.mbtiles \
    --use-attribute-for-id=feature_id \
    -zg \
    --minimum-zoom=4 \
    --read-parallel \
    --drop-densest-as-needed \
    --extend-zooms-if-still-dropping \
    --simplification=10 \
    --simplify-only-low-zooms \
    landslide_forest.geojson
```

NB that `--read-parallel` works with [GeoJSONSeq](https://gdal.org/en/latest/drivers/vector/geojsonseq.html)
(line-delimited GeoJSON with one feature per line and no wrapping FeatureCollection).
You could use `ogr2ogr` or something like this Python script to convert from regular
GeoJSON to a line-delimited series of features:

```python
import json
with open('features.geojson', 'r') as fh_in:
    with open('features.geojsonld', 'w') as fh_out:
        for f in data['features']:
            line = json.dumps(f)
            fh_out.write(line)
            fh_out.write("\n")
```

Once generated, the mbtiles file needs to sit in `./tileserver/vector/data` and have
an entry in [config.json or config-dev.json](https://github.com/nismod/infra-risk-vis/blob/c95b7cdcacf784fd92353f10f5f1312cfd0a5b6b/containers/vector/config.json#L11-L13).

The file and volume mapping for vector tiles is configured in the docker-compose
(e.g. [here for dev](https://github.com/nismod/infra-risk-vis/blob/c95b7cdcacf784fd92353f10f5f1312cfd0a5b6b/docker-compose-dev.yaml#L69)).

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

To bring up the database, refer to the [readme](../README.md) in the parent
directory for a full explanation of the required env files, etc., but briefly:

```bash
docker compose -f docker-compose-dev.yaml up db -d
```

### Awkward files

Unfortunately, not all the source data is openly available on the internet. Some
files must be manually copied into the appropriate location prior to running the
ETL process.

The affected datasets include:

- iris

Check the `snakemake` rules for more information, but source raster data should
typically reside in `raster/raw/<dataset>/`.

You may wish to remove write permissions to these files once they have been
installed, e.g. `chmod ug-w raster/raw/iris/*.tif`. This means `rm -r
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

### Running other rules by name

The top-level `Snakefile` imports rules from the pipelines using `module`
declarations, e.g.:

```
module land_cover:
    snakefile: "pipelines/land_cover/rules.smk"
    config: config
use rule * from land_cover as land_cover_*
```

To run any rule defined in a pipeline by name, include
the prefix specified in the line `use rule * from xxx as xxx_*`.

For example, the rule `download_300m_2020_from_CDS` from
`pipelines/land_cover/rules.smk` is imported `as land_cover_*`
and can be run using:

```bash
snakemake --cores 1 land_cover_download_300m_2020_from_CDS
```

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

- `"domain"` must be a short `lower_snake_case` string that will be exposed
  through the API and used by clients to request tiles. It will be prefixed by
  `terracotta_` to give the metadata database name. A good default choice would
  be the same string as the directory name for the pipeline, which is picked up
  as the `DATASET` wildcard by snakemake.
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

#### Removing a source from the tileserver metastore

Connect to the database server and delete the metadata database and reference row in `raster_tile_sources`:

```bash
# Drop a metadata database
psql -h localhost -U global_dev -c 'DROP DATABASE terracotta_storm;'
# Delete a row from raster tile sources
psql -h localhost -U global_dev -c "DELETE FROM raster_tile_sources WHERE domain = 'storm';"
```

Check the current state of the local database server:

```bash
# List all databases
psql -h localhost -U global_dev -c '\l'
# List raster tile sources
psql -h localhost -U global_dev -c 'select id, domain, name, "group", keys FROM raster_tile_sources;'
```

### Vector

To be implemented and then documented!
