import pandas as pd


def url_from_key(wildcards):
    """
    Lookup an Aqueduct TIFF URL from our layers file by KEY wildcard.
    """
    df: pd.DataFrame = pd.read_csv("pipelines/aqueduct/layers.csv")
    layer = df[df.filename == f"{wildcards.KEY}.tif"].squeeze()
    return layer.url


rule download_raw_data:
    """
    Download files from remote location.
    """
    input:
        layers = "pipelines/aqueduct/layers.csv",
    params:
        url = url_from_key
    output:
        "raster/raw/aqueduct/{KEY}.tif"
    shell:
        """
        wget {params.url} --output-document={output}
        """
