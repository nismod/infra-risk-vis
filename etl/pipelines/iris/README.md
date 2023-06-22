# Hazard - IRIS

Source notes - sby Nathan Sparks and Ralf Toumi, Imperial

## Pipeline

Data files stored locally within `/ouce-home/projects/mistral/iris`

-   Source Tiff Directories
-   Generate File Metadata from filenames & data dictionary
-   Snakemake
    -   extract from NetCDF to COG
    -   set zeros to nodata
    -   Move to output
-   Ingest to Terracotta DB

### API Metadata

```json
{
    "source_db": "iris",
    "global_type": "Hazard",
    "domain": "cyclone",
    "full_name": "Hazard Tropical Storm",
    "description": "description",
    "license": "license",
    "variables": {}
}
```
