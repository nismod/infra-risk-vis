import re

import pandas as pd
import xarray as xr


def netcdf_path_from_tiff(wildcards) -> pd.Series:
    """
    Lookup an IRIS source netCDF from our layers file by raster key.
    """
    fname = f"{wildcards.FILENAME}.tif"
    df: pd.DataFrame = pd.read_csv("pipelines/iris/layers.csv")
    layer = df[df.filename == fname].squeeze()
    return f"raster/raw/iris/{layer.nc_path}"


rule extract_netcdf_to_tiff:
    """
    Extract a return-period band from source IRIS netCDF files.
    """
    input:
        netcdf = netcdf_path_from_tiff
    output:
        tiff = "raster/raw/iris/{FILENAME}.tif"
    run:
        # extract return period (sub-string) from key wildcard
        rp, = re.search(r"rp_(\d+)", wildcards.FILENAME).groups()
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
