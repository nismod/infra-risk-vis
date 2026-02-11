set -e
set -x

dataset=$1

python scripts/dataset_tiffs_to_zarr.py \
    --dataset $dataset \
    --metadata-path pipelines/$dataset/metadata.json \
    --layers-path pipelines/$dataset/layers.csv \
    --layers-dir raster/cog/$dataset \
    --output-path raster/stack.zarr \
    && touch raster/stack/$dataset.flag
