# NbS Adaptation options

From within the `infra-risk-vis/etl` directory.

## Inputs from open-gira

Boundaries, HydroBASINS (with GADM codes) and the NbS opportunity areas are
open-gira outputs. Read them from a published open-gira results run (see
`etl/README.md`, "Vector layers from open-gira results"):

```bash
snakemake -s pipelines/nbs-adaptation/Snakefile -c1 --config \
    'open_gira_results={"base": "s3://nismod/irv-datapkg", "run": "20261124-fc04e93-r1"}'
```

or set `open_gira_results` in `pipelines/nbs-adaptation/config.yml`. With no
run configured, the same rules read a copy of an open-gira results directory
at `legacy_results_dir` (default
`/mnt/linux-filestore/mert2014/projects/open-gira/results`), replacing the
manual copies this README used to describe:

```bash
snakemake -s pipelines/nbs-adaptation/Snakefile -c1 \
    --config legacy_results_dir=/path/to/open-gira/results
```

| Output | From results layer | Legacy source |
| --- | --- | --- |
| `raw_data/boundaries/admin-level-{0,1,2}.geoparquet` | `admin_boundaries` (`admin_level`) | `input/admin-boundaries/admin-level-{0,1,2}.geoparquet` |
| `raw_data/boundaries/hybas_lev12_v1c_with_gadm_codes.geoparquet` | `hydrobasins_gadm` (`level=12`) | `input/hydrobasins/hybas_lev12_v1c_with_gadm_codes.geoparquet` |
| `raw_data/nbs-adaptation/{nbs}_with_EAD.geoparquet` | `nbs_opportunities` (`nbs`, `osm_extract`) | `slices/{osm_extract}_nbs/slice-*/{nbs}_with_EAD.parquet`, concatenated |

## Boundaries ETL

The Snakefile replaces `extract-nbs-adaptation-boundaries.ipynb`
(`scripts/nbs_adaptation.py`): rules `boundary`, `boundary_geojsonseq` and
`boundary_mbtiles`. The last calls [`tippecanoe`](https://github.com/felt/tippecanoe),
which needs to be installed or compiled from source.

Outputs:

- `raw_data/boundaries/{layer}.geojsonld` multi/polygons with added bbox_wkt property
- `raw_data/boundaries/{layer}_points.geojsonld` points from centroids
- `{tiles_dir}/{adm0,adm1,adm2,hybas}.mbtiles` (default `tiles_dir`:
  `../tileserver/vector/data`)
- `{tiles_dir}/{adm0,adm1,adm2,hybas}_points.mbtiles`, built from the points
  file (the notebook linked both MBTiles to the polygon file)

## Features ETL

Run the Snakefile (rule `nbs_opportunities`) for
`raw_data/nbs-adaptation/{nbs}_with_EAD.geoparquet`, then the
`extract-nbs-adaptation-opportunities.ipynb` notebook from its "Mangrove"
section (skip its "Extract" cell, which read the slices from the cluster).

- the Snakefile collects all opportunity areas (joined with EAD) from the
  results run (or the open-gira analysis "slices"), zero-filling the hazard
  columns the notebook uses, and saves combined GeoParquet files
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
