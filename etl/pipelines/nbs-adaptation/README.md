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

Run the `extract-nbs-adaptation.ipynb` notebook.

This calls [`tippecanoe`](https://github.com/felt/tippecanoe) to build MBTiles.
It will need to be installed or compiled from source.

Outputs:

- `raw_data/boundaries/{$l}.geojsonld` multi/polygons with added bbox_wkt property
- `raw_data/boundaries/{$l}_points.geojsonld` points from centroids
- `tileserver/vector/{adm0,adm1,adm2,hybas}.mbtiles`

## Features ETL

```bash
# put features in place
ls data/raw/btn__ls_nbs_current__split_ead.geojson
# run script
R --no-restore --no-save < run.R
```

Outputs:

- `btn_features_distinct.geojson`
- gather all non-hazard columns apart from feature_id and geometry into json
  column properties, write to CSV for load to database
- gather hazard columns, cost and benefit as adaptation cost benefit

```bash
tippecanoe -zg -o data/out/nbs.mbtiles data/out/btn_features_distinct.geojson --force
```
