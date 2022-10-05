# Hazard - STORM

Source URL (Present): https://data.4tu.nl/articles/dataset/STORM_tropical_cyclone_wind_speed_return_periods/12705164/3

Source URl (Future): https://data.4tu.nl/articles/dataset/STORM_climate_change_tropical_cyclone_wind_speed_return_periods/14510817/3

Trello Link:  https://trello.com/c/HjLKJJRB

## Pipeline

Manually Download Tiffs from Source -> `present` and `future` folders

```
Source Tiff Directories -> Generate File Metadata from filenames & data dictionary -> Snakemake [zero_nodata -> clip_north_south -> cloud_optimise -> Move to output] -> Ingest to Terracotta DB
```
