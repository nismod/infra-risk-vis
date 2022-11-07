# Hazard - Aqueduct

Source URL: http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/index.html

## Pipeline

```
Step 1: Source URL List -> Scrape for Tiffs -> Generate File Metadata from filenames & data dictionary -> Download files

Step 2: Snakemake [zero_nodata -> clip_north_south -> cloud_optimise -> Move to output]

Step 3: Ingest to Terracotta DB
```

### API Metadata

```json
[
	{
		"source_db": "aqueduct",
		"global_type": "Hazard",
		"domain": "fluvial",
		"full_name": "Hazard Aqueduct - Fluvial",
		"description": "description",
		"license": "license",
		"variables": {}
	},
	{
		"source_db": "aqueduct",
		"global_type": "Hazard",
		"domain": "coastal",
		"full_name": "Hazard Aqueduct - Coastal",
		"description": "description",
		"license": "license",
		"variables": {
			"some": "vars"
		}
	}
]
```