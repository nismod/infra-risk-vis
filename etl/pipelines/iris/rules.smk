import re

import pandas as pd
import xarray as xr


configfile: "../../config.yml"


def netcdf_path_from_key(wildcards) -> pd.Series:
    """
    Lookup an IRIS source netCDF from our layers file by raster key.
    """
    df: pd.DataFrame = pd.read_csv("pipelines/iris/layers.csv")
    layer = df[df.key == wildcards.KEY].squeeze()
    return f"raster/raw/iris/{layer.path}"


rule extract_netcdf_to_tiff:
    """
    Extract a return-period band from source IRIS netCDF files.
    """
    input:
        netcdf = netcdf_path_from_key
    output:
        tiff = "raster/raw/iris/{KEY}.tif"
    run:
        # extract return period (sub-string) from key wildcard
        rp, = re.search(r"rp_(\d+)", wildcards.KEY).groups()
        with xr.open_dataset(input.netcdf) as ds:
            ds.coords["longitude"].attrs = {
                "standard_name": "longitude",
                "long_name": "longitude",
                "units": "degrees_east",
                "axis": "X",
            }
            ds.coords["latitude"].attrs = {
                "standard_name": "latitude",
                "long_name": "latitude",
                "units": "degrees_north",
                "axis": "Y",
            }
            ds.rio.write_crs("epsg:4326", inplace=True)
            ds["vmax"].sel(rp=int(rp)).rio.to_raster(output.tiff)


def all_cog_file_paths(wildcards):
    """
    Generate list of IRIS cyclone wind speed output file paths.
    """
    df: pd.DataFrame = pd.read_csv("pipelines/iris/layers.csv")
    return expand("raster/cog/iris/{key}.tif", key=df.key)


rule ingest_rasters:
    """
    Create a dataset table in the MySQL database and ingest the cloud-optimised
    rasters to Terracotta.

    Requires the `tiles-db` MySQL service to be running.
    """
    input:
        all_cog_file_paths,
        script = "scripts/ingest.py",
        layers = "pipelines/iris/layers.csv",
        db_field_to_csv_header_map = "pipelines/iris/db_field_to_csv_header_map.json",
        tile_keys = "pipelines/iris/tile_keys.json",
    output:
        flag = "pipelines/iris/ingested_to_mysql.flag"
    shell:
        """
        python {input.script} load_csv \
            --internal_raster_base_path raster/cog/iris \
            --input_csv_filepath {input.layers} \
            --csv_to_db_field_map_path {input.db_field_to_csv_header_map} \
            --tile_keys_path {input.tile_keys} \
            --database_name iris

        touch {output.flag}
        """