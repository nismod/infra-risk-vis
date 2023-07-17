import pandas as pd


rule download_motorised:
    """
    Download motorised travel time to nearest healthcare facility.
    """
    output:
        archive = temp("raster/raw/traveltime_to_healthcare/202001_Global_Motorized_Travel_Time_to_Healthcare_2019.zip"),
        raster = "raster/raw/traveltime_to_healthcare/202001_Global_Motorized_Travel_Time_to_Healthcare_2019.tif"
    shell:
        """
        curl \
            'https://data.malariaatlas.org/geoserver/Accessibility/ows?service=CSW&version=2.0.1&request=DirectDownload&ResourceId=Accessibility:202001_Global_Motorized_Travel_Time_to_Healthcare' \
            -H 'Accept-Encoding: gzip, deflate, br' \
            -H 'Connection: keep-alive' \
            --output {output.archive}

        unzip {output.archive} $(basename {output.raster}) -d $(dirname {output.raster})
        """


rule download_walking:
    """
    Download walking travel time to nearest healthcare facility.
    """
    output:
        archive = temp("raster/raw/traveltime_to_healthcare/202001_Global_Walking_Only_Travel_Time_To_Healthcare_2019.zip"),
        raster = "raster/raw/traveltime_to_healthcare/202001_Global_Walking_Only_Travel_Time_To_Healthcare_2019.tif"
    shell:
        """
        curl \
            'https://data.malariaatlas.org/geoserver/Accessibility/ows?service=CSW&version=2.0.1&request=DirectDownload&ResourceId=Accessibility:202001_Global_Walking_Only_Travel_Time_To_Healthcare' \
            -H 'Accept-Encoding: gzip, deflate, br' \
            -H 'Connection: keep-alive' \
            --output {output.archive}

        unzip {output.archive} $(basename {output.raster}) -d $(dirname {output.raster})
        """


def all_cog_file_paths(wildcards):
    """
    Generate list of travel time to healthcare facility output file paths.
    """
    df: pd.DataFrame = pd.read_csv("pipelines/traveltime_to_healthcare/layers.csv")
    return expand("raster/cog/traveltime_to_healthcare/{key}.tif", key=df.key)


rule ingest_rasters:
    """
    Create a dataset table in the MySQL database and ingest the cloud-optimised
    rasters to Terracotta.

    Requires the `tiles-db` MySQL service to be running.
    """
    input:
        all_cog_file_paths,
        script = "scripts/ingest.py",
        layers = "pipelines/traveltime_to_healthcare/layers.csv",
        db_field_to_csv_header_map = "pipelines/traveltime_to_healthcare/db_field_to_csv_header_map.json",
        tile_keys = "pipelines/traveltime_to_healthcare/tile_keys.json",
    output:
        flag = "pipelines/traveltime_to_healthcare/ingested_to_mysql.flag"
    shell:
        """
        python {input.script} load_csv \
            --internal_raster_base_path raster/cog/traveltime_to_healthcare \
            --input_csv_filepath {input.layers} \
            --csv_to_db_field_map_path {input.db_field_to_csv_header_map} \
            --tile_keys_path {input.tile_keys} \
            --database_name traveltime_to_healthcare

        touch {output.flag}
        """