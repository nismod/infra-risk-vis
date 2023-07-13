import pandas as pd


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


def all_cog_file_paths(wildcards):
    """
    Generate list of ESA land cover output file paths.
    """
    df: pd.DataFrame = pd.read_csv("pipelines/esa_land_cover/layers.csv")
    return expand("raster/cog/esa_land_cover/{key}.tif", key=df.key)


rule ingest_rasters:
    """
    Create a dataset table in the MySQL database and ingest the cloud-optimised
    rasters to Terracotta.

    Requires the `tiles-db` MySQL service to be running.
    """
    input:
        all_cog_file_paths,
        script = "scripts/ingest.py",
        layers = "pipelines/esa_land_cover/layers.csv",
        db_field_to_csv_header_map = "pipelines/esa_land_cover/db_field_to_csv_header_map.json",
        tile_keys = "pipelines/esa_land_cover/tile_keys.json",
    output:
        flag = "pipelines/esa_land_cover/ingested_to_mysql.flag"
    shell:
        """
        python {input.script} load_csv \
            --internal_raster_base_path raster/cog/esa_land_cover \
            --input_csv_filepath {input.layers} \
            --csv_to_db_field_map_path {input.db_field_to_csv_header_map} \
            --tile_keys_path {input.tile_keys} \
            --database_name esa_land_cover

        touch {output.flag}
        """

