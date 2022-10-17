# Hazard - Aqueduct

Source URL: http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/index.html

Trello Link:  https://trello.com/c/k7tVeD3E

## Pipeline

```
Step 1: Source URL List -> Scrape for Tiffs -> Generate File Metadata from filenames & data dictionary -> Download files

Step 2: Snakemake [zero_nodata -> clip_north_south -> cloud_optimise -> Move to output]

Step 3: Ingest to Terracotta DB
```