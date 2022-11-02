# JRC World Population

Sourced from here:  https://ghsl.jrc.ec.europa.eu/download.php?ds=pop

## Pipeline

All Snakefile:

tif -> reproject -> clip north south -> COG

### API Metadata

```json
{
  "source_db": "jrc_pop",
  "global_type": "Exposure",
  "domain": "population",
  "full_name": "JRC Population",
  "description": "description",
  "license": "license",
  "variables": {}
}
```
