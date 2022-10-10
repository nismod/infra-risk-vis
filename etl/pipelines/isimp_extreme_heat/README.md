# Hazard - ISIMP Extreme Heat

Source URL: https://data.isimip.org/datasets/4f79e7aa-7def-4665-854f-93ff033bec37/

Trello Link:  https://trello.com/c/hqKMGs55

## Pipeline

World Pop 2020 Sourced from: https://data.humdata.org/dataset/worldpop-population-counts-for-world/resource/677d30ab-896e-44e5-9a31-05452bc3124b

Resampled to 0.5deg using:  gdalwarp -t_srs EPSG:4326 -tr 0.5 0.5 -r near -te_srs EPSG:4326 -of GTiff pipelines/isimp_extreme_heat/ppp_2020_1km_Aggregated.tif pipelines/isimp_extreme_heat/05deg.tiff