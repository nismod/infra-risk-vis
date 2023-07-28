rule download_300m_2020_from_CDS:
    """
    Download ESA land cover classification data from Copernicus CDS.

    Subsequently extract and process to TIFF.
    """
    output:
        archive = temp("raster/raw/esa_land_cover/download.zip"),
        netcdf = temp("raster/raw/esa_land_cover/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc"),
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


rule generate_terracotta_colourmap:
    input:
        csv = "pipelines/esa_land_cover/colourmap.csv",
        script = "scripts/legend_to_tc_colourmap.py"
    output:
        json = "pipelines/esa_land_cover/colourmap.json"
    shell:
        """
        python {input.script} {input.csv} {output.json}
        """


rule ingest_categorical_raster:
    """
    Custom ingestion rule for ESA land cover as raster is categorical and needs
    a mapping from raster integer value to colour and classification.
    """
    input:
        raster = "raster/cog/esa_land_cover/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.tif",
        script = "scripts/ingest.py",
        csv_colourmap = "pipelines/esa_land_cover/colourmap.csv",
        tile_key_mapping = "pipelines/esa_land_cover/tile_key_mapping.json",
    output:
        flag = "pipelines/esa_land_cover/rasters_ingested.flag"
    shell:
        """
        python {input.script} load_single_categorical \
            --categorical_legend_csv_filepath {input.csv_colourmap} \
            --categorical_csv_label_column Label \
            --categorical_csv_value_column Value \
            --database_name esa_land_cover \
            --input_raster_filepath {input.raster} \
            --internal_raster_base_path raster/cog/esa_land_cover \
            --categorical_key_values_json_path {input.tile_key_mapping}

        touch {output.flag}
        """