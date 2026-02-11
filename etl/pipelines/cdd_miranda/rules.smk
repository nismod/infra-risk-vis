"""Changes in Cooling Degree Days (CDD) between the 1.5ºC and 2.0ºC global warming scenarios
"""
rule download:
    output:
        nc_abs="raster/raw/cdd_miranda/CDD_absolute_mean_change_from_15_to_20.nc",
        nc_rel="raster/raw/cdd_miranda/CDD_relative_mean_change_from_15_to_20.nc",
        nc_std="raster/raw/cdd_miranda/CDD_stdv_change_from_15_to_20.nc",
    shell:
        """
        wget -nc https://ora.ox.ac.uk/objects/uuid:8d95c423-816c-4a4f-88b6-eb7a040cb40e/files/d8910jv156 \
            --output-document=raster/raw/cdd_miranda/CDD_absolute_mean_change_from_15_to_20.nc
        wget -nc https://ora.ox.ac.uk/objects/uuid:8d95c423-816c-4a4f-88b6-eb7a040cb40e/files/dj38607487 \
            --output-document=raster/raw/cdd_miranda/CDD_relative_mean_change_from_15_to_20.nc
        wget -nc https://ora.ox.ac.uk/objects/uuid:8d95c423-816c-4a4f-88b6-eb7a040cb40e/files/d6w924c429 \
            --output-document=raster/raw/cdd_miranda/CDD_stdv_change_from_15_to_20.nc
        """

rule extract:
    input:
        netcdf="raster/raw/cdd_miranda/CDD_{ABS_REL}_mean_change_from_15_to_20.nc",
    output:
        tiff="raster/raw/cdd_miranda/CDD_{ABS_REL}_mean_change_from_15_to_20.tif"
    run:
        import xarray as xr
        with xr.open_dataset(input.netcdf) as ds:
            ds.rio.write_crs("epsg:4326", inplace=True)
            ds.coords["longitude0"].attrs = {
                "standard_name": "longitude",
                "long_name": "longitude",
                "units": "degrees_east",
                "axis": "X",
            }
            ds.coords["latitude0"].attrs = {
                "standard_name": "latitude",
                "long_name": "latitude",
                "units": "degrees_north",
                "axis": "Y",
            }
            ds["CDD_total"].rio.to_raster(output.tiff)
