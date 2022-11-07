# Hazard - ISIMP Drought

Trello Link:  https://trello.com/c/r7gJNfSV

## Pipeline

JRC Populaton sourced from: https://ghsl.jrc.ec.europa.eu/download.php?ds=pop

1. Python Script isimp_drought.py 

__NOTE__: This matches the extreme heat script - they could / shout be merged into a more general ISIMP data handler

Download List File -> Generate meta -> generate occurrence (temporally averaged) tiffs -> generate exposure (world pop * avg) tiffs -> generate hazard_csv

2. Snakemake

Zero nodata -> clip north/south -> export to Cloud Optimised GeoTiff

### API Metadata

```json
{
  "source_db": "drought",
  "global_type": "Hazard",
  "domain": "drought",
  "full_name": "ISIMP Drought",
  "description": "description",
  "license": "license",
  "variables": {}
}
```