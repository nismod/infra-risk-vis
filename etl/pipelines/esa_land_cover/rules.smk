rule download_300m_2020_from_CDS:
    """
    Download ESA land cover classification data from Copernicus CDS.

    Subsequently extract and process to TIFF.
    """
    output:
        archive = temp("raster/raw/esa_land_cover/download.zip"),
        netcdf = "raster/raw/esa_land_cover/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc",
        tif = "raster/raw/esa_land_cover/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.tif",
    resources:
        mem_mb = 10000,
        disk_mb = 4000
    shell:
        """
        python pipelines/esa_land_cover/download_from_CDS.py \
            --dataset_name satellite-land-cover \
            --variable all \
            --file_format zip \
            --version v2.1.1 \
            --year 2020 \
            --output_path {output.archive}

        unzip {output.archive} $(basename {output.netcdf}) -d $(dirname {output.netcdf})

        gdalwarp \
            -of Gtiff \
            -co COMPRESS=LZW \
            -ot Byte \
            -te -180.0000000 -90.0000000 180.0000000 90.0000000 \
            -tr 0.002777777777778 0.002777777777778 \
            -t_srs EPSG:4326 \
            NETCDF:{output.netcdf}:lccs_class \
            {output.tif}
        """