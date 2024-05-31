"""These rules define interactions with the terracotta metadata databases that
store raster / tileset state.

N.B. Dataset specific rules are located in pipelines/<dataset>.
"""

import pandas as pd


rule ingest_rasters:
    """Create a dataset table in the terracotta metadata database and ingest the
    cloud-optimised rasters.
    """
    input:
        tiffs=lambda wildcards: expand(
            "raster/cog/{dataset}/{filename}",
            dataset=wildcards.DATASET,
            filename=pd.read_csv(f"pipelines/{wildcards.DATASET}/layers.csv").filename
        ),
        layers = "pipelines/{DATASET}/layers.csv",
        metadata = "pipelines/{DATASET}/metadata.json",
    output:
        flag = "raster/ingest/{DATASET}.flag"
    script:
        "scripts/ingest.py"


rule ingest_metadata:
    """Creates a row in the main database
    """
    input:
        ingest_flag = "raster/ingest/{DATASET}.flag",
        metadata = "pipelines/{DATASET}/metadata.json",
    output:
        flag = "raster/metadata/{DATASET}.flag"
    script:
        "scripts/ingest_metadata.py"
