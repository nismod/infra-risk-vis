# Extract-transform-load

## Requirements
- Postgres database
  - schema defined in `../backend/backend/db/models.py`
- Python, snakemake and other packages
  - requirements in `../backend/Pipfile`
- [`jq`](https://stedolan.github.io/jq/)
- [`geojson-polygon-labels`](https://github.com/andrewharvey/geojson-polygon-labels)
- [`tippecanoe`](https://github.com/mapbox/tippecanoe)

### GDAL tools

[ogr2ogr](https://www.gdal.org/ogr2ogr.html) and other GDAL programs are used
for spatial data processing. On Ubuntu, run:

    sudo apt-get install gdal-bin

### Tippecanoe

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
