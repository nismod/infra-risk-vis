# ESA Landcover

Categorical data originally supplied in netcdf format: `C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc`

__NOTE__: This dataset fails to load into a Cloud-hosted MySQL database due to its size (the connection in Terracotta drops before the ingest completes.)  The work-around is to load it to a local database, then replicate the database remotely using `mysqldump`.

## Pipeline

Extract lccs_class from netcdf -> clip -> COG

### API Metadata

```json
{
  "source_db": "land_cover",
  "global_type": "Exposure",
  "domain": "land_cover",
  "full_name": "ESA Land Cover",
  "description": "description",
  "license": "license",
  "variables": {}
}
```

