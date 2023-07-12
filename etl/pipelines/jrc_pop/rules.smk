import pandas as pd


configfile: "../../config.yml"


def url_from_key(wildcards):
    """
    Lookup a JRC population URL from our metadata file by KEY wildcard.
    """
    df: pd.DataFrame = pd.read_csv("pipelines/jrc_pop/metadata.csv")
    metadata = df[df.key == wildcards.KEY].squeeze()
    return metadata.url


rule download_and_unzip_raw_data:
    """
    Download JRC population data from remote location and unzip it.
    """
    input:
        "pipelines/jrc_pop/metadata.csv"
    params:
        url = url_from_key
    output:
        zip_file = temp("raster/raw/jrc_pop/{KEY}.zip"),
        raster = "raster/raw/jrc_pop/{KEY}.tif",
    resources:
        disk_mb=1000
    shell:
        """
        wget {params.url} --output-document={output.zip_file}
        unzip {output.zip_file} $(basename {output.raster}) -d $(dirname {output.raster})
        """


def gdalwarp_bounds(bbox: dict[str, float]) -> str:
    """
    Given dict of `minx`, `miny`, `maxx` and `maxy`, return coordinates as
    concatentated string in correct order for gdalwarp -te argument.
    """
    extents = (bbox["minx"], bbox["miny"], bbox["maxx"], bbox["maxy"])
    return " ".join(map(lambda extent: f"{extent:.3f}", extents))


rule clip_and_reproject_raster:
    """
    Reproject from Mollweide to WGS84 and clip to bounds.
    """
    input:
        raster = "raster/no_data/jrc_pop/{KEY}.tif",
        target_CRS = "pipelines/jrc_pop/WGS84_CRS.wkt",
    params:
        bounds = gdalwarp_bounds(config["raster_bounds"])
    output:
        "raster/clip/jrc_pop/{KEY}.tif"
    resources:
        mem_mb=10000
    priority:
        80,
    shell:
        """
        gdalwarp \
            -t_srs {input.target_CRS} \
            -of GTiff \
            -co COMPRESS=LZW \
            -te {params.bounds} \
            {input.raster} \
            {output}
        """


def all_cog_file_paths(wildcards):
    """
    Generate list of JRC population output file paths.
    """
    df: pd.DataFrame = pd.read_csv("pipelines/jrc_pop/metadata.csv")
    return expand("raster/cog/jrc_pop/{key}.tif", key=df.key)


rule ingest_rasters:
    """
    Create a dataset table in the MySQL database and ingest the cloud-optimised
    rasters to Terracotta.

    Requires the `tiles-db` MySQL service to be running.
    """
    input:
        all_cog_file_paths,
        script = "scripts/ingest.py",
        metadata = "pipelines/jrc_pop/metadata.csv",
        db_field_to_csv_header_map = "pipelines/jrc_pop/db_field_to_csv_header_map.json",
        tile_keys = "pipelines/jrc_pop/tile_keys.json",
    output:
        flag = "pipelines/jrc_pop/ingested_to_mysql.flag"
    shell:
        """
        python {input.script} load_csv \
            --internal_raster_base_path raster/cog/jrc_pop \
            --input_csv_filepath {input.metadata} \
            --csv_to_db_field_map_path {input.db_field_to_csv_header_map} \
            --tile_keys_path {input.tile_keys} \
            --database_name jrc_pop

        touch {output.flag}
        """
