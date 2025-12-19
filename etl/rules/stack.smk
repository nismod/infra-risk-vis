"""Create a combined Zarr store containing groups for each dataset
"""

import pandas as pd


rule raster_dataset_tiffs_to_zarr:
    """Create a group in the combined Zarr store with the tensor
    of n-dimensional raster data (typically 2D latitude/longitude,
    plus epoch, climate scenario, and other model variables).
    """
    input:
        script = "scripts/dataset_tiffs_to_zarr.py",
        tiffs=lambda wildcards: expand(
            "{tiff_dir}/{filename}",
            tiff_dir=f"raster/cog/{wildcards.DATASET}",
            filename=pd.read_csv(f"pipelines/{wildcards.DATASET}/layers.csv").filename
        ),
        layers = "pipelines/{DATASET}/layers.csv",
        metadata = "pipelines/{DATASET}/metadata.json",
    params:
        zarr = "raster/stack.zarr",
        tiff_dir = "raster/cog/{DATASET}",
    output:
        flag = "raster/stack/{DATASET}.flag",
    shell:
        """
        python {input.script} \
            --dataset {wildcards.DATASET} \
            --metadata-path {input.metadata} \
            --layers-path {input.layers} \
            --layers-dir {params.tiff_dir} \
            --output-path {params.zarr} \
            && touch {output.flag}
        """
