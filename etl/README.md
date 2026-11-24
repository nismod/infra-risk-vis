# Extract Transform Load (ETL) for infra-risk-vis (IRV)

## Introduction

This directory contains the logic, configuration and documentation for
acquiring, processing and ingesting data such that a running IRV stack may host
as tilesets.

The data processing steps are broadly as follows:

### Raster

- Download data if possible, or fetch it from the irv-datapkg catalogue (see
  [Sourcing data from irv-datapkg](#sourcing-data-from-irv-datapkg))
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

## Sourcing data from irv-datapkg

[irv-datapkg](https://github.com/nismod/irv-datapkg) is becoming the single
producer of analysis-ready (T1) global open data, published as a static STAC
catalogue under one base URI and pinned by six-monthly releases. Once a
dataset is published there, its pipeline here stops downloading from upstream
and fetches the catalogue file instead.

The cut point is the raw raster path. `fetch_datapkg` (`rules/datapkg.smk`)
writes exactly `raster/raw/{pipeline}/{key}.tif` for every row of the
pipeline's `layers.csv`, verified against its sha256. Everything from
`clip_raster` onwards (Web Mercator clipping, cloud optimisation, terracotta,
Zarr stack, database) is unchanged and stays in this repository (irv-datapkg
ADR 0007). T1 rasters keep the source CRS and grid (ADR 0001), so
`clip_raster`, which already warps to EPSG:4326, does any reprojection
(e.g. GHSL Mollweide, SoilGrids Homolosine).

### Switching a pipeline over

Switching is per pipeline and off by default. In `config.yml` (or with
`--config`):

```yaml
datapkg:
  base: s3://example-bucket/irv-datapkg   # catalogue base URI
  private_base: null                      # default: {base}-private
  release: "2026.1"
  lock: datapkg.lock.yaml                 # used instead of base/release once it exists
  cache_dir: null                         # default: $IRV_DATAPKG_CACHE; neither: no cache
  enabled: [aqueduct, storm]
```

For each enabled pipeline, `fetch_datapkg` claims only the raw paths listed in
its `layers.csv` and takes precedence (`ruleorder`) over every pipeline rule
writing to `raster/raw/{pipeline}/`. Other pipelines, and other raw paths, keep
their own download rules, which remain as the fallback until each dataset is
published. Check what will run with:

```bash
snakemake -n raster/metadata/aqueduct.flag
```

The client streams each file straight to `raster/raw/`, verifying size and
sha256 while streaming and renaming into place atomically (no second copy on
disk). Set `cache_dir` (or `IRV_DATAPKG_CACHE`) to keep a download cache
instead: raw files are then hardlinked from it and are read-only.

The fetch runs as a shell rule, so credentials for a private or S3 catalogue
can go in `../envs/{environment}/.etl.env` (e.g. `AWS_PROFILE`, or
`IRV_DATAPKG_CACHE`). The client needs `irv-datapkg>=0.3.0` (in
`environment.yml`; 0.3.0 is not tagged yet, so it installs from git); with
nothing enabled the workflow does not import it.

### Mapping layers to catalogue items

Each pipeline that can be switched has a `pipelines/{pipeline}/datapkg.yml`:

```yaml
dataset: aqueduct_flood          # catalogue dataset id
properties:                      # layers.csv column -> STAC item property
  hazard: hazard
  rp: return_period
  rcp: rcp
  epoch: epoch
  gcm: gcm
values:                          # layers.csv value -> property value
  hazard: {fluvial: river}
  gcm: {None: null}              # null: property absent or null
filters: {}                      # extra filters for every row
filename: "{id}.tif"             # for rows generated from new items
columns: {var: "depth_{hazard}"} # other columns of generated rows
```

Each `layers.csv` row resolves to exactly one item, by the optional
`datapkg_item` column (item id) if set, otherwise by the property filter built
from the mapping and the optional `datapkg_filter` column
(`key=value;key=v1|v2`). A `null` value (in `datapkg_filter`, or a `values`
map to `null`) matches a property that is absent or null, as `None` /
`Catalog.ABSENT` does in the client's `Catalog.items()`; resolution uses
`items()` itself, so from a lock file items whose recorded properties rule
them out are not even read. An optional `datapkg_dataset`
column overrides `dataset` for pipelines that combine several catalogue
datasets (e.g. `nature`). Filenames (the raster keys used by terracotta and
the frontend) never change.

Recording `datapkg_item` is recommended: it makes resolution exact; with a
lock file the item is fetched from the lock entry alone (`fetch_locked`: href,
sha256 and size, no STAC reads), only the named items are locked, and the
item's checksum becomes a job parameter, so a new release re-runs only the
layers whose bytes changed. Fill it from the catalogue with:

```bash
python scripts/layers_from_catalogue.py aqueduct            # refresh: fill datapkg_item
python scripts/layers_from_catalogue.py aqueduct --add-missing  # also add rows for new items
python scripts/layers_from_catalogue.py aqueduct --check    # CI: non-zero if stale/unmatched
```

`--add-missing` (or a pipeline with no `layers.csv`) generates rows from item
properties, reversing the `values` maps and rendering `filename`/`columns`.
Catalogue options default to the `datapkg` section of `config.yml`;
`--base/--release/--lock` override them.

### Locking to a release

Consumers lock to a release, never "latest" (irv-datapkg ADR 0003):

```bash
snakemake --cores 1 datapkg_lock   # writes datapkg.lock.yaml for the enabled pipelines
git add datapkg.lock.yaml
```

`datapkg_lock` runs `irv-datapkg lock --base ... --release ...` with
`--item DATASET/ID` for every layer that names a `datapkg_item`, so only the
items this deployment serves are locked; a dataset in which any layer is still
resolved by filter is locked whole. While the lock file exists it is used
instead of `release`, and the client refuses any item whose checksum differs
from it. To move to a new release, change `release`, delete the lock and
re-run `datapkg_lock`.

In CI, check that the committed lock is current (same release and versions,
nothing retired, every needed item locked with its current checksum) without
writing anything:

```bash
snakemake --cores 1 datapkg_lock_check
# equivalent to:
irv-datapkg lock --check --base "$BASE" --release 2026.1 \
    --output datapkg.lock.yaml --item aqueduct_flood/ITEM_ID ... [DATASET ...]
```

#### Mirrors

The lock records the catalogue `base` it was made from. When `datapkg.base` is
also set, it deliberately overrides the lock's base (`fetch_datapkg` passes
both `--base` and `--lock`; `Catalog.from_lock(lock, base=...)`): the locked
versions, item ids and checksums still come from the lock, but files are read
from that base, e.g. a mirror or a local copy of the catalogue. The private
base is then `datapkg.private_base` or `{base}-private`. Leave `base` null, or
equal to the recorded base, to read from where the lock was made.

### metadata.json from the STAC collection

```bash
python scripts/metadata_from_catalogue.py aqueduct          # update license + description
python scripts/metadata_from_catalogue.py aqueduct --check  # CI
```

This sets `license` (SPDX id) and `description` from the collection.
`raster_tile_sources` has no citation column and the frontend shows
`description` as the source credit, so `description` takes `sci:citation`
(plus DOI); `--description description` uses the collection description
instead. `domain`, `name`, `group`, `keys` and colour maps stay local.

### Pipeline mapping status

The `datapkg.yml` files are provisional until each dataset is published: the
property names follow the catalogue spec, and values must be confirmed against
the published items (then run `layers_from_catalogue.py`).

| Pipeline | Catalogue dataset(s) | Status | layers.csv column -> property | Replaces / notes |
| --- | --- | --- | --- | --- |
| aqueduct | `aqueduct_flood` | wired | hazard (fluvial->river), rp->return_period, rcp, epoch, gcm (None->null) | `download_raw_data`. Coastal items will add subsidence / sea-level percentile dimensions: add to `filters` or `datapkg_filter` |
| storm | `storm_fixed_rp` | wired | rp->return_period, gcm, rcp, epoch | `download` (Zenodo zip) + `unpack` |
| ghsl_pop | `ghsl_pop` | wired | epoch; filter resolution=1000 | `download_and_unzip_raw_data`. Mollweide grid kept at T1 |
| jrc_flood | `jrc_flood_glofas` | wired | rp->return_period | Assumes one global mosaic COG per RP from irv-datapkg; the local rules never produced these files |
| landslide_arup | `arup_landslide` | wired | subtype->(dimension TBD by WP6) | `download_metadata` + `download_landslides` |
| ghsl_buildings | `ghsl_built_s` | todo | subtype (all/non_residential) -> e.g. `class`; filter epoch=2020, resolution=1000 | `download_all_building_types`, `download_non_residential_type` |
| gem_earthquake | `gem_seismic_hazard` | todo | rp->return_period, medium->site | `download` + `unpack`. Single item |
| isimip | `isimip_heat_drought` | todo | hazard, metric, rcp, epoch, gcm, impact_model | `download` + `unpack`. Upstream version still to reconcile (two Zenodo records) |
| land_cover | `copernicus_lulc` | todo | land_cover -> item id (single layer) | `download_300m_2020_from_CDS` (CDS API, netCDF->GeoTIFF moves to irv-datapkg). Keep categorical ingest |
| dem | `merit_dem` | todo | derivative->variable (elevation/slope) | `download` |
| cdd_miranda | `cdd_miranda` | todo | metric (absolute/relative) | `download` + `extract` (netCDF->COG moves to irv-datapkg) |
| nature | `soilgrids_soc`, `biodiversity_intactness`, `forest_landscape_integrity` | todo | `datapkg_dataset` per row; single item each | Three download rules. Homolosine kept at T1 (`clip_raster` warps). Confirm the FLII 0-10000 -> 0-10 rescale is applied at T1 or via `raster:bands` scale |
| social | `relative_wealth_index`, `hdi_mosaiks` | todo | `datapkg_dataset` per row | hdi: `download_gridded_hdi`. rwi: T1 is points (GeoParquet), so `process_rwi_grid` rasterisation stays here, reading the fetched parquet rather than a raw tif |
| traveltime_to_healthcare | `traveltime_healthcare` | todo | subtype->mode (motorized/walking) | `download_motorised`, `download_walking` |
| iris | `iris` (private) | todo | epoch, rp->return_period, ssp | Manual copy + `extract_netcdf_to_tiff`. Needs `private_base` credentials |
| osm_rail, osm_roads, gem_exposure, nbs-adaptation | `openstreetmap`, ... | out of scope | n/a | Vector / separate Snakefiles, not `raster/raw` |

### Retiring the download rules

The pipeline download rules listed under "Replaces / notes" above are
deprecated but kept as the fallback for pipelines that are not enabled. A
pipeline's download rules (and `unpack`/`extract` steps) are removed only when:

1. its catalogue dataset version is pinned by a published irv-datapkg release,
   and the remote passes `irv-datapkg validate --deep` and `irv-datapkg audit
   --release`;
2. the pipeline is in the committed `datapkg.enabled`, `datapkg.lock.yaml`
   pins that release, and `snakemake --cores 1 datapkg_lock_check` passes;
3. a run on the production machine with the pipeline enabled reproduced the
   raster COGs and terracotta metadata, recorded in the removal PR;
4. any note in the table (dimensions to confirm, rescaling, rasterisation) is
   resolved.

The cross-repository list of legacy rules and their status is kept in
irv-datapkg [`docs/deprecation.md`](https://github.com/nismod/irv-datapkg/blob/main/docs/deprecation.md).
`ci/etl-datapkg.yml` (to copy to `.github/workflows/`) runs the tests below in CI.

### Tests

```bash
pip install snakemake pytest pandas xarray -e path/to/irv-datapkg
cd etl && python -m pytest tests
```

The tests build a tiny synthetic catalogue with the irv-datapkg producer API
(nothing is downloaded), check the mapping and layers generation, dry-run the
workflow with pipelines enabled and disabled, and run `datapkg_lock` and
`fetch_datapkg` for real in a scratch working directory.

See also [docs/serving-bundle.md](docs/serving-bundle.md) for publishing the
processed serving bundle.

## Adding new datasets

### Raster

To add additional raster datasets to the ETL pipeline you will need several new
files in `pipelines/<new_dataset_name>`:

- `README.md` describing the dataset
- `layers.csv` table of layers
- `metadata.json` metadata to be written to the Postgres database
- `rules.smk` to acquire and process the data into a WGS84 raster. The rule(s)
  should write output files to `raster/raw/<dataset>/`. If the dataset is
  published in irv-datapkg, write a `datapkg.yml` instead (see
  [Sourcing data from irv-datapkg](#sourcing-data-from-irv-datapkg)) and keep
  any download rules only as a fallback. The shared rules will
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

#### Deploying a raster layer

After testing against a local backend, `snakemake` can run the final ingest/metadata
steps against the production database.

Remove the flag files, then run with `--config environment=prod` to source environment
files for the production database.

```
rm raster/{ingest,metadata}/social.flag
snakemake --cores 1 --config environment=prod  -- raster/metadata/social.flag
```

Copy up the COGs and update the `backend` service (see the top-level [README](../README.md)).

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
