rule download:
    """
    Fetch our preprocessed data from zenodo.

    N.B. Requesting any file will trigger a download of the whole archive
    """
    output:
        "raster/raw/isimip/{KEY}.tif",
    shell:
        """
        OUTPUT_DIR=$(dirname {output})
        ARCHIVE=$OUTPUT_DIR/lange2020_expected_occurrence.zip

        zenodo_get --record=8147088 --output-dir=$OUTPUT_DIR

        unzip $ARCHIVE -d $OUTPUT_DIR

        mv $OUTPUT_DIR/data/* $OUTPUT_DIR
        rm -r $OUTPUT_DIR/data
        """


rule POST_metadata_to_backend:
    """
    Requires the `backend` and postgreSQL `db` services to be running.

    This rule is a special case as we have two objects to create, one for
    extreme heat and one for drought. They both share a database on the MySQL
    instance, but require their own metadata store in postgreSQL.
    """
    input:
        ingest_flag = "raster/ingest/isimip.flag",
        heat_metadata = "pipelines/isimip/metadata_extreme_heat.json",
        drought_metadata = "pipelines/isimip/metadata_drought.json",
    output:
        flag = "raster/metadata/isimip.flag"
    shell:
        """
        http --check-status --follow POST http://$GATEWAY_HOST:$GATEWAY_PORT/api/tiles/sources x-token:$BE_API_TOKEN < {input.heat_metadata}
        http --check-status --follow POST http://$GATEWAY_HOST:$GATEWAY_PORT/api/tiles/sources x-token:$BE_API_TOKEN < {input.drought_metadata}

        touch {output.flag}
        """
