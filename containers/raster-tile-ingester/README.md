# Raster Tileserver Ingester

Utility for managing Cloud-Optimised Tiff files in a Terracotta MySQL database.

This script can work with a single categorical raster, or a CSV containing multiple rasters to be loaded (non-categorical).

### Example CSV Load

__NOTE__: This can be achieved using `docker-compose-dev.yaml` -> `raster-tile-ingester` utility container.

1.  Prepare a CSV containing information about the rasters to be loaded:

```csv
hazard,metric,path,rcp,epoch,gcm,key
drought,occurrence,/mnt/d/oxford/isimp_drought/occurrence/lange2020_clm45_gfdl-esm2m_ewembi_rcp26_2005soc_co2_led_global_annual_2006_2099_2030_occurrence.tif,2x6,2030,gfdl-esm2m,lange2020_clm45_gfdl-esm2m_ewembi_rcp26_2005soc_co2_led_global_annual_2006_2099_2030_occurrence
drought,occurrence,/mnt/d/oxford/isimp_drought/occurrence/lange2020_clm45_gfdl-esm2m_ewembi_rcp26_2005soc_co2_led_global_annual_2006_2099_2050_occurrence.tif,2x6,2050,gfdl-esm2m,lange2020_clm45_gfdl-esm2m_ewembi_rcp26_2005soc_co2_led_global_annual_2006_2099_2050_occurrence
drought,occurrence,/mnt/d/oxford/isimp_drought/occurrence/lange2020_clm45_gfdl-esm2m_ewembi_rcp26_2005soc_co2_led_global_annual_2006_2099_2080_occurrence.tif,2x6,2080,gfdl-esm2m,lange2020_clm45_gfdl-esm2m_ewembi_rcp26_2005soc_co2_led_global_annual_2006_2099_2080_occurrence
```

2. Prepare the ingest script arguments

```bash
ingest.py load_csv \ # The load operation
--internal_raster_base_path /data/traveltime_to_healthcare # The top-level internal path under-which raster entries will be stored in the DB (e.g. a local path of /home/me/raster.tiff will be stored as /data/traveltime_to_healthcare/raster.tiff)
--input_csv_filepath /opt/hazard_layers.csv # Path to CSV containing rasters and keys (see above)
--tile_keys type,rcp,epoch,gcm # ORDERED list of tile keys to be used - these must be included as headers in the CSV file.
--csv_key_column_map "{\"file_basename\": \"key\", \"type\": \"hazard\", \"rcp\": \"rcp\", \"epoch\": \"epoch\", \"gcm\": \"gcm\"}" # Naming of columns in the DB can be changed from what is included in the CSV by using this map.
--database_name traveltime_to_healthcare # The database to be used (will be created if not already present)
```

## Usage

```bash
usage: ingest.py [-h] [--input_csv_filepath INPUT_CSV_FILEPATH] [--input_raster_filepath INPUT_RASTER_FILEPATH]
                 [--categorical_legend_csv_filepath CATEGORICAL_LEGEND_CSV_FILEPATH]
                 [--categorical_csv_label_column CATEGORICAL_CSV_LABEL_COLUMN] [--categorical_csv_value_column CATEGORICAL_CSV_VALUE_COLUMN]
                 [--categorical_key_values CATEGORICAL_KEY_VALUES] [--tile_keys TILE_KEYS] [--csv_key_column_map CSV_KEY_COLUMN_MAP]
                 [--database_name DATABASE_NAME] [--internal_raster_base_path INTERNAL_RASTER_BASE_PATH]
                 {load_csv,load_single_categorical,delete_database_entries,drop_database}

Terracotta Ingester

positional arguments:
  {load_csv,load_single_categorical,delete_database_entries,drop_database}
                        Type of load operation - - `load_csv` Loads all rasters in a CSV file generated from an ETL pipeline, e.g. with the
                        header: `hazard,metric,path,rcp,epoch,gcm,key` - `load_single_categorical` Loads a single categorical raster (using
                        provided category_map). Requires: `input_raster_filepath`, `categorical_legend_csv_filepath`,
                        `categorical_csv_label_column`, `categorical_csv_value_column` - `delete_database_entries` Delete all raster entries
                        from a given database (leaving the database empty). Requires `database_name` - `drop_database` Drop a database and all
                        raster entries contained-within. Requires `database_name`

options:
  -h, --help            show this help message and exit
  --input_csv_filepath INPUT_CSV_FILEPATH
                        Absolute path to the CSV file containing information for each raster
  --input_raster_filepath INPUT_RASTER_FILEPATH
                        Absolute path to a categorical raster file being loaded
  --categorical_legend_csv_filepath CATEGORICAL_LEGEND_CSV_FILEPATH
                        Absolute path to the CSV file containing legend information about a categorical raster being loaded
  --categorical_csv_label_column CATEGORICAL_CSV_LABEL_COLUMN
                        Name of the label column in CSV
  --categorical_csv_value_column CATEGORICAL_CSV_VALUE_COLUMN
                        Name of the value column in CSV
  --categorical_key_values CATEGORICAL_KEY_VALUES
                        Valid JSON string containing key:value mapping for the categorical raster. NOTE: this will be loaded as an OrderedDict
                        - so the key ordering will be maintained for DB usage. e.g.
                        {"type":"exposure","sensor":"S2","data":"20181010","band":"cloudmask"}
  --tile_keys TILE_KEYS
                        A comma-seperated list of tile keys and ordering, e.g. hazard,rp,rcp,gcm
  --csv_key_column_map CSV_KEY_COLUMN_MAP
                        Map of DB Keys to column names for input CSV (must be valid JSON str and contain key: file_basename), e.g.
                        '{"file_basename": "key", "type": "hazard", "rp": "rp", "rcp": "rcp", "epoch": "epoch", "gcm": "gcm"}'
  --database_name DATABASE_NAME
                        Name of the output database (will be created if it doesnt exist)
  --internal_raster_base_path INTERNAL_RASTER_BASE_PATH
                        Internal path to the raster file
```


## Environment

The following environment variables are required:

```
#### Terracotta Env
TC_DRIVER_PATH=mysql://{mysql user}:{mysql passwd}@{mysql host}
TC_DRIVER_PROVIDER=mysql
TC_PNG_COMPRESS_LEVEL=0
TC_RESAMPLING_METHOD="nearest"
TC_REPROJECTION_METHOD="nearest"
```

## Docker

A utility container for loading data can be build from the included `Dockerfile`
