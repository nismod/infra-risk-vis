import re

import pandas as pd
import xarray as xr


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