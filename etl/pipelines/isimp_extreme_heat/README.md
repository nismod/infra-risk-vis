# Hazard - ISIMP Extreme Heat

Source URL: https://data.isimip.org/datasets/4f79e7aa-7def-4665-854f-93ff033bec37/

Trello Link:  https://trello.com/c/hqKMGs55

## Pipeline

JRC Populaton sourced from: https://ghsl.jrc.ec.europa.eu/download.php?ds=pop

[ Python Script isimp_extreme_heat.py ]

Download List File -> Generate meta -> generate occurrence (temporally averaged) tiffs -> generate exposure (world pop * avg) tiffs -> generate hazard_csv

[ Snakemake ]

Zero nodata -> clip north/south -> export to Cloud Optimised GeoTiff
