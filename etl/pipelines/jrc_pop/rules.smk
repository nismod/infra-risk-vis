import pandas as pd

from pipelines.helpers import gdalwarp_bounds

configfile: "../../config.yml"


def url_from_key(wildcards):
    """
    Lookup a JRC population URL from our layers file by KEY wildcard.
    """
    df: pd.DataFrame = pd.read_csv("pipelines/jrc_pop/layers.csv")
    layer = df[df.key == wildcards.KEY].squeeze()
    return layer.url


rule download_and_unzip_raw_data:
    """
    Download JRC population data from remote location and unzip it.
    """
    input:
        "pipelines/jrc_pop/layers.csv"
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


rule clip_and_reproject_raster:
    """
    Reproject from Mollweide to WGS84 and clip to bounds.
    """
    input:
        raster = "raster/no_data/jrc_pop/{KEY}.tif",
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
            -t_srs EPSG:4326 \
            -of GTiff \
            -co COMPRESS=LZW \
            -te {params.bounds} \
            {input.raster} \
            {output}
        """