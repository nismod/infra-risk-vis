# Loading OSM Roads and Rail Output from OpenGIRA into Infra-Risk Database

Geoparquet -> PG -> geojsonl -> MBTiles.

Features are kept in PG

Feature type is added to PG FeatureLayers table

ID in Features table is auto-increment

string_id in features table is the OSM id (not osm_way_id - which is duplicated)

## Prerequisites:

- Geoparquet OSM files including ways
- Addition of required feature-classes (one per row) in `network_tilelayers.csv` (`asset_type` is used to down-select on parquet load - which maps to `assert_type_column` in `network_layers.csv`)
- Mapping of requires feature-classes to their source files in `network_layers.csv`
- PG Environment variables in the run-environment (e.g. in snakemake .env file if using docker-compose)

No other pre-processing is required