import pandas as pd


def all_cog_file_paths(wildcards):
    """
    Generate list of GEM earthquake output file paths.
    """
    df: pd.DataFrame = pd.read_csv("pipelines/gem_earthquake/layers.csv")
    return expand("raster/cog/gem_earthquake/{key}.tif", key=df.key)


rule ingest_rasters:
    """
    Create a dataset table in the MySQL database and ingest the cloud-optimised
    rasters to Terracotta.

    Requires the `tiles-db` MySQL service to be running.
    """
    input:
        all_cog_file_paths,
        script = "scripts/ingest.py",
        layers = "pipelines/gem_earthquake/layers.csv",
        db_field_to_csv_header_map = "pipelines/gem_earthquake/db_field_to_csv_header_map.json",
        tile_keys = "pipelines/gem_earthquake/tile_keys.json",
    output:
        flag = "pipelines/gem_earthquake/ingested_to_mysql.flag"
    shell:
        """
        python {input.script} load_csv \
            --internal_raster_base_path raster/cog/gem_earthquake \
            --input_csv_filepath {input.layers} \
            --csv_to_db_field_map_path {input.db_field_to_csv_header_map} \
            --tile_keys_path {input.tile_keys} \
            --database_name gem_earthquake

        touch {output.flag}
        """