# GEM Exposure Vector to MBTiles

GEM data came in two packages:

1. Geopackage of Seismic Exposure per country

__NOTE__: This gpkg file was missing a number of features globally - including Morocco, Arunachal Predesh and Kashmir

2. CSV of flooding exposure per country

The above datasets were joined on `WB_A3` (the only common column apart from `WB_NAME`).  

__NOTE__: The geopackage contained 11 additional features which did not link (nulls for flooding exposure)

## Pipeline

Geopackage -> Joined to Flooding eposure CSV -> Exported as `geojsonl` -> Snakemake to `mbtile` (via tippecanoe)
