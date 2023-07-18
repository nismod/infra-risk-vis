import pandas as pd


configfile: "../../config.yml"


checkpoint create_hazard_csv_file:
    """
    Create table of Aqueduct hazard layers.

    N.B. We use `shell` rather than `script` to allow for command line
    arguments.  However, using `shell` with the script path hard coded would not
    trigger a re-run on the modification of the script. Instead treat the script
    as an input to retain this behaviour.
    """
    input:
        script = "pipelines/aqueduct/parser.py"
    params:
        source_url = "http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2",
        list_url = "http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/index.html"
    output:
        csv = "pipelines/aqueduct/layers.csv"
    shell:
        """
        python {input.script} --list_url {params.list_url} --source_url {params.source_url} --csv_path {output.csv}
        """


def url_from_key(wildcards):
    """
    Lookup an Aqueduct TIFF URL from our layers file by KEY wildcard.
    """
    df: pd.DataFrame = pd.read_csv(checkpoints.create_hazard_csv_file.get(**wildcards).output.csv)
    layer = df[df.key == wildcards.KEY].squeeze()
    return layer.url


rule download_raw_data:
    """
    Download files from remote location.
    """
    input:
        layers = lambda wildcards: checkpoints.create_hazard_csv_file.get().output,
    params:
        url = url_from_key
    output:
        "raster/raw/aqueduct/{KEY}.tif"
    shell:
        """
        wget {params.url} --output-document={output}
        """


def all_cog_file_paths(wildcards):
    """
    Generate list of Aqueduct output file paths.
    """
    df: pd.DataFrame = pd.read_csv(checkpoints.create_hazard_csv_file.get(**wildcards).output.csv)
    return expand("raster/cog/aqueduct/{key}.tif", key=df.key)


rule ingest_rasters:
    """
    Create a dataset table in the MySQL database and ingest the cloud-optimised
    rasters to Terracotta.

    Requires the `tiles-db` MySQL service to be running.

    TODO: Can this rule be factored out into main Snakefile? Currently, it
    depends too tightly on the create_hazard_csv_file checkpoint to do so.
    """
    input:
        all_cog_file_paths,
        script = "scripts/ingest.py",
        layers = "pipelines/aqueduct/layers.csv",
        db_field_to_csv_header_map = "pipelines/aqueduct/db_field_to_csv_header_map.json",
        tile_keys = "pipelines/aqueduct/tile_keys.json",
    output:
        flag = "pipelines/aqueduct/ingested_to_mysql.flag"
    shell:
        """
        python {input.script} load_csv \
            --internal_raster_base_path raster/cog/aqueduct \
            --input_csv_filepath {input.layers} \
            --csv_to_db_field_map_path {input.db_field_to_csv_header_map} \
            --tile_keys_path {input.tile_keys} \
            --database_name aqueduct

        touch {output.flag}
        """


rule POST_metadata_to_backend:
    """
    Requires the `backend` and postgreSQL `db` services to be running.

    This rule is a special case as we have two objects to create, one for fluvial and one for coastal.
    """
    input:
        ingest_flag = "pipelines/aqueduct/ingested_to_mysql.flag",
        fluvial_metadata = "pipelines/aqueduct/metadata_fluvial.json",
        coastal_metadata = "pipelines/aqueduct/metadata_coastal.json",
    output:
        flag = "pipelines/aqueduct/posted_to_backend.flag"
    shell:
        """
        # N.B. 4XX responses result in a zero-valued httpie exit status
        http POST http://$BE_HOST:$BE_PORT/tiles/sources x-token:$BE_API_TOKEN < {input.fluvial_metadata}
        http POST http://$BE_HOST:$BE_PORT/tiles/sources x-token:$BE_API_TOKEN < {input.coastal_metadata}

        touch {output.flag}
        """