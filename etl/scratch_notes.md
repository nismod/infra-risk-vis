# Data updates notes

Some things done outside of the etl pipeline, or rerun.


Drop damages and adaptation results before reloading for irrigation.
```sql
delete from damages_rp where feature_id >= 400000000 and feature_id < 600000000;
delete from damages_expected where feature_id >= 400000000 and feature_id < 600000000;
delete from damages_npv where feature_id >= 400000000 and feature_id < 600000000;
delete from adaptation_cost_benefit where feature_id >= 400000000 and feature_id < 600000000;
```

Convert to GeoJSON
```bash
ogr2ogr -f GeoJSONSeq marine_combined.geojsonl marine_combined.gpkg
```

Attempt using snakemake defaults to run mbtiles
```bash
snakemake -c20 --configfile=config_otherdata.yml "../tileserver/vector/data/drought_options.mbtiles"
```

Custom values for NbS tile layer generation
```bash
tippecanoe \
    --use-attribute-for-id=uid \
    --read-parallel \
    --output=../tileserver/vector/data/drought_combined.mbtiles \
    --layer=drought_combined \
    --drop-smallest-as-needed \
    -zg \
    --force \
    vector/drought_combined.geojsonl

tippecanoe \
    --use-attribute-for-id=uid \
    --read-parallel \
    --output=../tileserver/vector/data/natural_marine_combined.mbtiles \
    --layer=natural_marine_combined \
    --drop-smallest-as-needed \
    -zg \
    --force \
    vector/natural_marine_combined.geojsonl

tippecanoe \
    -t tmp \
    --generate-ids \
    --read-parallel \
    --output=../tileserver/vector/data/natural_terrestrial_combined.mbtiles \
    --layer=natural_terrestrial_combined \
    --drop-smallest-as-needed \
    --extend-zooms-if-still-dropping \
    -zg \
    --force \
    vector/natural_terrestrial_combined.geojsonl
```

Run split/convert/combine for terrestrial as points tilelayer (too big for
ogr2ogr in one go)
```bash
split -n l/16 vector/natural_terrestrial_combined.geojsonl

ls x* | parallel -j20 'ogr2ogr \
    -sql \
    "SELECT ST_Centroid(geometry), * FROM {}" \
    -dialect sqlite \
    {}_points.geojsonl \
    {}'

cat x*geojsonl > vector/natural_terrestrial_combined_points.geojsonl
# use -t tmp to use relative ./tmp folder for temp storage if local disk /tmp is full
tippecanoe \
    -t tmp \
    --generate-ids \
    --read-parallel \
    --output=../tileserver/vector/data/natural_terrestrial_combined_points.mbtiles \
    --layer=natural_terrestrial_combined_points \
    --drop-densest-as-needed \
    --extend-zooms-if-still-dropping \
    -zg \
    --force \
    vector/natural_terrestrial_combined_points.geojsonl
```
