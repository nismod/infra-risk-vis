# NbS Adaptation options

From within the `infra-risk-vis/etl` directory:

```bash
# incoming boundaries and NbS opportunities from open-gira
mkdir -p raw_data/boundaries
mkdir -p raw_data/nbs-adaptation
# load-ready data as MBTiles for vector tileserver and CSVs for database
mkdir -p vector/csvs

# Copy prepared/extracted boundaries from open-gira
export OG_RESULTS_DIR="/mnt/linux-filestore/mert2014/projects/open-gira/results"
cp $OG_RESULTS_DIR/input/admin-boundaries/admin-level-0.geoparquet ./raw_data/boundaries
cp $OG_RESULTS_DIR/input/admin-boundaries/admin-level-1.geoparquet ./raw_data/boundaries
cp $OG_RESULTS_DIR/input/admin-boundaries/admin-level-2.geoparquet ./raw_data/boundaries
cp $OG_RESULTS_DIR/input/hydrobasins/hybas_lev12_v1c_with_gadm_codes.geoparquet ./raw_data/boundaries
```

## Boundaries ETL

Run the `extract-nbs-adaptation-boundaries.ipynb` notebook.

This calls [`tippecanoe`](https://github.com/felt/tippecanoe) to build MBTiles.
It will need to be installed or compiled from source.

Outputs:

- `raw_data/boundaries/{$l}.geojsonld` multi/polygons with added bbox_wkt property
- `raw_data/boundaries/{$l}_points.geojsonld` points from centroids
- `tileserver/vector/data/{adm0,adm1,adm2,hybas}.mbtiles`
- `tileserver/vector/data/{adm0,adm1,adm2,hybas}_points.mbtiles`

## Features ETL

Run the `extract-nbs-adaptation-opportunities.ipynb` notebook.

- collects all opportunity areas (joined with EAD) from open-gira analysis
  "slices", saves to combined GeoParquet files
- gathers all non-hazard columns apart from feature_id and geometry into json
  column properties, write to CSV for load to database
- gathers hazard columns, cost and benefit as adaptation cost benefit

Outputs:

- GeoParquet and GPKG for archive/data analysis:
  - `raw_data/nbs-adaptation/landslide_slope_vegetation_with_EAD{.geoparquet,_grouped.geoparquet,_grouped_gt0.gpkg}`
  - `raw_data/nbs-adaptation/mangrove_with_EAD{.geoparquet,_grouped.geoparquet,_grouped_gt0.gpkg}`
  - `raw_data/nbs-adaptation/river_basin_afforestation_with_EAD{.geoparquet,_grouped.geoparquet,_grouped_gt0.gpkg}`
- CSVs for load to database:
  - `raw_data/nbs-adaptation/{cf,ls,rf}_{features,adaptation_cost_benefit}.csv`
- GeoJSON for transform to MBTiles:
  - `raw_data/nbs-adaptation/nbs_{cf,ls,rf}.json`
  - `raw_data/nbs-adaptation/nbs_{cf,ls,rf}_points.json`
- MBTiles:
  - `tileserver/vector/data/nbs_{cf,ls,rf}.mbtiles`
  - `tileserver/vector/data/nbs_{cf,ls,rf}_points.mbtiles`
