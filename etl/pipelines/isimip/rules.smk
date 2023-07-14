rule download_isimip:
    """
    Fetch our preprocessed data from zenodo.
    """
    output:
        # requesting any file will trigger a download of the whole archive
        requested_raster = "raster/raw/isimip/{KEY}.tif",
    shell:
        """
        zenodo_get --record=8147088 --output-dir=$(dirname {output.requested_raster})
        """


def all_cog_file_paths(wildcards):
    """
    Generate list of ISIMIP extreme heat / drought output file paths.
    """
    df: pd.DataFrame = pd.read_csv("pipelines/isimip/layers.csv")
    return expand("raster/cog/isimip/{key}.tif", key=df.key)


rule ingest_rasters:
    """
    Create a dataset table in the MySQL database and ingest the cloud-optimised
    rasters to Terracotta.

    Requires the `tiles-db` MySQL service to be running.
    """
    input:
        all_cog_file_paths,
        script = "scripts/ingest.py",
        layers = "pipelines/isimip/layers.csv",
        db_field_to_csv_header_map = "pipelines/isimip/db_field_to_csv_header_map.json",
        tile_keys = "pipelines/isimip/tile_keys.json",
    output:
        flag = "pipelines/isimip/ingested_to_mysql.flag"
    shell:
        """
        python {input.script} load_csv \
            --internal_raster_base_path raster/cog/isimip \
            --input_csv_filepath {input.layers} \
            --csv_to_db_field_map_path {input.db_field_to_csv_header_map} \
            --tile_keys_path {input.tile_keys} \
            --database_name isimip

        touch {output.flag}
        """