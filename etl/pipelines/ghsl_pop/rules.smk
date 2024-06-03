import pandas as pd


configfile: "../../config.yml"


def url_from_key(wildcards):
    """
    Lookup a JRC population URL from our layers file by KEY wildcard.
    """
    df: pd.DataFrame = pd.read_csv("pipelines/ghsl_pop/layers.csv")
    layer = df[df.filename == f"{wildcards.KEY}.tif"].squeeze()
    return layer.url


rule download_and_unzip_raw_data:
    """
    Download JRC population data from remote location and unzip it.
    """
    input:
        "pipelines/ghsl_pop/layers.csv"
    params:
        url = url_from_key
    output:
        zip_file = temp("raster/raw/ghsl_pop/{KEY}.zip"),
        raster = "raster/raw/ghsl_pop/{KEY}.tif",
    resources:
        disk_mb=1000
    shell:
        """
        wget {params.url} --output-document={output.zip_file}
        unzip {output.zip_file} $(basename {output.raster}) -d $(dirname {output.raster})
        """
