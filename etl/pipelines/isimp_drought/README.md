# Hazard - ISIMP Drought

Trello Link:  https://trello.com/c/r7gJNfSV

## Pipeline

JRC Populaton sourced from: https://ghsl.jrc.ec.europa.eu/download.php?ds=pop

[ Python Script isimp_drought.py ]

__NOTE__: This matches the extreme heat script - they could / shout be merged into a more general ISIMP data handler

Download List File -> Generate meta -> generate occurrence (temporally averaged) tiffs -> generate exposure (world pop * avg) tiffs -> generate hazard_csv

[ Snakemake ]

Zero nodata -> clip north/south -> export to Cloud Optimised GeoTiff
