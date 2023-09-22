"""
These rules define interactions with the MySQL and postgreSQL databases that
store raster / tileset state.

N.B. Dataset specific rules are located in pipelines/<dataset>.
"""

import pandas as pd


rule ingest_rasters:
    """
    Create a dataset table in the MySQL database and ingest the cloud-optimised
    rasters to Terracotta.

    Requires the `tiles-db` MySQL service to be running.
    """
    input:
        lambda wildcards: expand(
            "raster/cog/{dataset}/{key}.tif",
            dataset=wildcards.DATASET,
            key=pd.read_csv(f"pipelines/{wildcards.DATASET}/layers.csv").key
        ),
        script = "scripts/ingest.py",
        layers = "pipelines/{DATASET}/layers.csv",
        db_field_to_csv_header_map = "pipelines/{DATASET}/db_field_to_csv_header_map.json",
        tile_keys = "pipelines/{DATASET}/tile_keys.json",
    output:
        flag = "raster/ingest/{DATASET}.flag"
    shell:
        """
        python {input.script} load_csv \
            --local_raster_base_path raster/cog/{wildcards.DATASET} \
            --db_raster_base_path /data/{wildcards.DATASET} \
            --input_csv_filepath {input.layers} \
            --csv_to_db_field_map_path {input.db_field_to_csv_header_map} \
            --tile_keys_path {input.tile_keys} \
            --database_name {wildcards.DATASET}

        touch {output.flag}
        """


rule POST_metadata_to_backend:
    """
    Requires the `backend` and postgreSQL `db` services to be running.
    """
    input:
        ingest_flag = "raster/ingest/{DATASET}.flag",
        metadata = "pipelines/{DATASET}/metadata.json",
    output:
        flag = "raster/metadata/{DATASET}.flag"
    shell:
        """
        http --check-status --follow POST http://$GATEWAY_HOST:$GATEWAY_PORT/api/tiles/sources x-token:$BE_API_TOKEN < {input.metadata}

        touch {output.flag}
        """
