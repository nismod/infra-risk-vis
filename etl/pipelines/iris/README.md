# Hazard - IRIS

Authored by Nathan Sparks and Ralf Toumi, Imperial College London.

## Pipeline

Data files stored locally within `/ouce-home/projects/mistral/iris`.

Expect source netCDF files to be within `raster/raw/iris/`, unfortunately not yet publically available.

- Generate file metadata from filenames & data dictionary (included as layers.csv)
- Snakemake
    - Extract tiffs from NetCDF
    - Set zeros to nodata
    - Clip extent
    - Cloud optimise
    - Ingest to Terracotta MySQL DB
    - Create table metadata in PostgreSQL DB
