# Extract-transform-load

ETL is conducted using a combination of pre-processing scripts and Snakemake configurations to generate raster and vector tiles for hosting.

These pipelines broken-out into subdirectories under `./pipelines`

## Pipelines

[ISIMP Aqueduct (Coastal and Riverine Flooding)](pipelines/aqueduct/README.md)

[ESA Land Cover (Categorical Data on Global Land Cover)](pipelines/esa_land_cover/README.md)

[Natural Asset Exposure (Organic Carbon, Forest Integrity and Biodiversity Intactness)](pipelines/exposure_nature/README.md)

[GEM Earthquake  (Siesmic hazard data from GEM)](pipelines/gem_earthquake/README.md)

[GEM Exposure (Seismic and Flooding Exposure data from GEM)](pipelines/gem_exposure/README.md)

[GHSL Buildings (Single raster containing information about residential and non-residential buildings globally)](pipelines/ghsl_buildings/README.md)

[ISIMP Drought (ISIMP Drought hazard)](pipelines/isimp_drought/README.md)

[ISIMP Extreme Heat (ISIMP Extreme Heat)](pipelines/isimp_extreme_heat/README.md)

[JRC Population (Single raster containing information about global population)](pipelines/jrc_pop/README.md)

[OSM Rail (Open Streetmap vector data on rail edges and stations)](pipelines/osm_rail/README.md)

[OSM Roads (Open Streetmap vector data on road edges)](pipelines/osm_roads/README.md)

[STORM (Global cyclone rasters)](pipelines/storm/README.md)

[Traveltime to Healthcare (Shows amount of time to travel to healthcare globally)](pipelines/traveltime_to_healthcare/README.md)

## Setup - Docker

A Docker image can be built and used for running Snakemake workflows using `Dockerfile`.

This file also gives information about installing the required dependencies locally.

## Setup - Local Install

A large number of dependencies are required - see Dockerfile for more information.

### Requirements
- Postgres database
  - schema defined in `../backend/backend/db/models.py`
- Python, snakemake and other packages
  - requirements in `../backend/Pipfile`
- [`jq`](https://stedolan.github.io/jq/)
- [`geojson-polygon-labels`](https://github.com/andrewharvey/geojson-polygon-labels)
- [`tippecanoe`](https://github.com/mapbox/tippecanoe)

#### GDAL tools

[ogr2ogr](https://www.gdal.org/ogr2ogr.html) and other GDAL programs are used
for spatial data processing. On Ubuntu, run:

    sudo apt-get install gdal-bin

#### Tippecanoe

The data preparation steps use
[Mapbox tippecanoe](https://github.com/mapbox/tippecanoe) to build vector tiles
from large feature sets.

The easiest way to install tippecanoe on OSX is with Homebrew:

    brew install tippecanoe

On Ubuntu it will usually be easiest to build from the source repository:

    sudo apt-get install build-essential g++ libsqlite3-dev zlib1g-dev
    git clone https://github.com/mapbox/tippecanoe
    cd tippecanoe
    make -j
    make
