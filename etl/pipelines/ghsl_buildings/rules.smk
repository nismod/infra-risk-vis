import pandas as pd

from pipelines.helpers import gdalwarp_bounds

configfile: "../../config.yml"


rule download_all_building_types:
    """
    Download JRC built up area layer.
    """
    output:
        archive = temp("raster/raw/ghsl_buildings/all_types_E{EPOCH}_R{RELEASE}_archive.zip"),
        raster = "raster/raw/ghsl_buildings/GHS_BUILT_S_E{EPOCH}_GLOBE_R{RELEASE}_54009_1000_V1_0.tif"
    shell:
        """
        wget \
            https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R{wildcards.RELEASE}/GHS_BUILT_S_E{wildcards.EPOCH}_GLOBE_R{wildcards.RELEASE}_54009_1000/V1-0/GHS_BUILT_S_E{wildcards.EPOCH}_GLOBE_R{wildcards.RELEASE}_54009_1000_V1_0.zip \
            --output-document={output.archive}
        
        unzip {output.archive} $(basename {output.raster}) -d $(dirname {output.raster})
        """


rule download_non_residential_type:
    """
    Download JRC non-residential built up area layer.
    """
    output:
        archive = temp("raster/raw/ghsl_buildings/non_residential_E{EPOCH}_R{RELEASE}_archive.zip"),
        raster = "raster/raw/ghsl_buildings/GHS_BUILT_S_NRES_E{EPOCH}_GLOBE_R{RELEASE}_54009_1000_V1_0.tif"
    shell:
        """
        wget \
            https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R{wildcards.RELEASE}/GHS_BUILT_S_NRES_E{wildcards.EPOCH}_GLOBE_R{wildcards.RELEASE}_54009_1000/V1-0/GHS_BUILT_S_NRES_E{wildcards.EPOCH}_GLOBE_R{wildcards.RELEASE}_54009_1000_V1_0.zip \
            --output-document={output.archive}
        
        unzip {output.archive} $(basename {output.raster}) -d $(dirname {output.raster})
        """


rule clip_and_reproject_raster:
    """
    Reproject from Mollweide to WGS84 and clip to bounds.
    """
    input:
        raster = "raster/no_data/ghsl_buildings/{KEY}.tif",
    params:
        bounds = gdalwarp_bounds(config["raster_bounds"])
    output:
        "raster/clip/ghsl_buildings/{KEY}.tif"
    resources:
        mem_mb=10000
    priority:
        80,
    shell:
        """
        gdalwarp \
            -t_srs EPSG:4326 \
            -of GTiff \
            -co COMPRESS=LZW \
            -te {params.bounds} \
            {input.raster} \
            {output}
        """


def all_cog_file_paths(wildcards):
    """
    Generate list of GHSL built up settlement output file paths.
    """
    df: pd.DataFrame = pd.read_csv("pipelines/ghsl_buildings/layers.csv")
    return expand("raster/cog/ghsl_buildings/{key}.tif", key=df.key)


rule ingest_rasters:
    """
    Create a dataset table in the MySQL database and ingest the cloud-optimised
    rasters to Terracotta.

    Requires the `tiles-db` MySQL service to be running.
    """
    input:
        all_cog_file_paths,
        script = "scripts/ingest.py",
        layers = "pipelines/ghsl_buildings/layers.csv",
        db_field_to_csv_header_map = "pipelines/ghsl_buildings/db_field_to_csv_header_map.json",
        tile_keys = "pipelines/ghsl_buildings/tile_keys.json",
    output:
        flag = "pipelines/ghsl_buildings/ingested_to_mysql.flag"
    shell:
        """
        python {input.script} load_csv \
            --internal_raster_base_path raster/cog/ghsl_buildings \
            --input_csv_filepath {input.layers} \
            --csv_to_db_field_map_path {input.db_field_to_csv_header_map} \
            --tile_keys_path {input.tile_keys} \
            --database_name ghsl_buildings

        touch {output.flag}
        """