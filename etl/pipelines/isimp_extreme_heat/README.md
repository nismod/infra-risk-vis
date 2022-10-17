# Hazard - ISIMP Extreme Heat

Source URL: https://data.isimip.org/datasets/4f79e7aa-7def-4665-854f-93ff033bec37/

Trello Link:  https://trello.com/c/hqKMGs55

## Pipeline

World Pop 2020 Sourced from: https://data.humdata.org/dataset/worldpop-population-counts-for-world/resource/677d30ab-896e-44e5-9a31-05452bc3124b

JRC Populaton sourced from: https://ghsl.jrc.ec.europa.eu/download.php?ds=pop

Resampled to 0.5deg using:  gdalwarp -t_srs EPSG:4326 -te -180 -90 180 90 -tr 0.5 0.5 -r average -te_srs EPSG:4326 -of GTiff pipelines/isimp_extreme_heat/GHS_POP_P2030_GLOBE_R2022A_54009_1000_V1_0.tif pipelines/isimp_extreme_heat/jrc_popn_05deg.tiff