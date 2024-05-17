rule download:
    """
    Fetch our preprocessed data from zenodo.
    """
    output:
        zip="raster/raw/isimip/lange2020_expected_occurrence.zip",
    shell:
        """
        pushd $(dirname {output.zip})
            zenodo_get -w links.txt --record=8147088
            wget -nc -i links.txt
            md5sum -c md5sums.txt
        popd
        """

rule unpack:
    input:
        zip=rules.download.output.zip,
    output:
        tiff="raster/raw/isimip/{KEY}.tif",
    shell:
        """
        unzip -n {input.zip} $(basename {output.tiff}) -d $(dirname {output.tiff})
        """


rule POST_metadata_to_backend:
    """
    Requires the `backend` and postgreSQL `db` services to be running.

    This rule is a special case as we have two objects to create, one for
    extreme heat and one for drought. They both share a terracotta metadata database
    but require their own metadata store in postgreSQL.
    """
    input:
        ingest_flag = "raster/ingest/isimip.flag",
        heat_metadata = "pipelines/isimip/metadata_extreme_heat.json",
        drought_metadata = "pipelines/isimip/metadata_drought.json",
    output:
        flag = touch("raster/metadata/isimip.flag")
    shell:
        """
        curl -X POST \
            -H 'Content-Type: application/json' \
            -H 'X-Token: $BE_API_TOKEN' \
            -d @{input.heat_metadata} \
            http://$GATEWAY_HOST:$GATEWAY_PORT/api/tiles/sources

        curl -X POST \
            -H 'Content-Type: application/json' \
            -H 'X-Token: $BE_API_TOKEN' \
            -d @{input.drought_metadata} \
            http://$GATEWAY_HOST:$GATEWAY_PORT/api/tiles/sources

        touch raster/metadata/isimip.flag
        """
