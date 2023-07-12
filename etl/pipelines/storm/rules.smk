def all_cog_file_paths(wildcards):
    """
    Generate list of STORM tropical cyclone output file paths.
    """
    df: pd.DataFrame = pd.read_csv("pipelines/storm/layers.csv")
    return expand("raster/cog/storm/{key}.tif", key=df.key)


rule ingest_rasters:
    """
    Create a dataset table in the MySQL database and ingest the cloud-optimised
    rasters to Terracotta.

    Requires the `tiles-db` MySQL service to be running.
    """
    input:
        all_cog_file_paths,
        script = "scripts/ingest.py",
        layers = "pipelines/storm/layers.csv",
        db_field_to_csv_header_map = "pipelines/storm/db_field_to_csv_header_map.json",
        tile_keys = "pipelines/storm/tile_keys.json",
    output:
        flag = "pipelines/storm/ingested_to_mysql.flag"
    shell:
        """
        python {input.script} load_csv \
            --internal_raster_base_path raster/cog/storm \
            --input_csv_filepath {input.layers} \
            --csv_to_db_field_map_path {input.db_field_to_csv_header_map} \
            --tile_keys_path {input.tile_keys} \
            --database_name storm

        touch {output.flag}
        """