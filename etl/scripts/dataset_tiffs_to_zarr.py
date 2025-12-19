"""
Docstring for etl.scripts.dataset_tiffs_to_zarr
"""

import json
import logging
import warnings
from pathlib import Path

import click
import dask
import pandas
import xarray as xr
from snail.io import extend_rasters_metadata
from zarr.errors import UnstableSpecificationWarning

warnings.filterwarnings(
    "ignore",
    message="Consolidated metadata is currently not part in the Zarr format*",
    category=UserWarning,
)
warnings.filterwarnings(
    "ignore",
    message="The data type*",
    category=UnstableSpecificationWarning,
)


@click.command()
@click.version_option("1.0.0")
@click.option(
    "--dataset",
    "-d",
    required=True,
    type=str,
    help="Dataset short name",
)
@click.option(
    "--metadata-path",
    "-m",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to dataset metadata JSON",
)
@click.option(
    "--layers-path",
    "-l",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to layers metadata CSV",
)
@click.option(
    "--layers-dir",
    "-ld",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="Directory to read layers TIFFS.",
)
@click.option(
    "--output-path",
    "-o",
    required=True,
    type=click.Path(
        exists=False, file_okay=False, dir_okay=True, readable=True, writable=True
    ),
    help="Directory to write Zarr store.",
)
def main(dataset, metadata_path, layers_path, layers_dir, output_path) -> None:
    """Rewrite TIFFs to Zarr store"""
    ll_chunksize = 2048
    logging.info("Converting dataset TIFFs to Zarr for %s", dataset)
    with open(metadata_path, "r") as fh:
        metadata = json.load(fh)
    layers = pandas.read_csv(layers_path)

    store_path = Path(output_path)
    if not store_path.exists():
        store_path.mkdir(mode=0o750, parents=False, exist_ok=True)

    source_path = Path(layers_dir)
    layer_paths = [source_path / fname for fname in layers.filename.tolist()]
    layers["path"] = layer_paths
    layers, _ = extend_rasters_metadata(layers)
    for col in metadata["keys"]:
        if layers[col].dtype == "object":
            layers.fillna({col: ""}, inplace=True)

    for var, grid_layers in layers.groupby("var"):
        grid_layer_paths = grid_layers.path.tolist()
        assert len(grid_layers.grid_id.unique()) == 1, f"{var=} had non-aligned grids"
        meta_dicts = grid_layers[metadata["keys"]].to_dict(orient="records")
        logging.info(
            f"Processing {len(grid_layer_paths)} layers for {dataset}/{var} in {source_path}"
        )
        setup_store(
            store_path,
            grid_layer_paths[0],
            dataset,
            var,
            grid_layers,
            metadata,
            ll_chunksize=ll_chunksize,
        )

        for layer_path, meta in zip(grid_layer_paths, meta_dicts):
            logging.info("Read %s", layer_path)
            layer_ds = (
                xr.open_dataset(
                    layer_path,
                    engine="rasterio",
                    chunks=dict(band=1, x=ll_chunksize, y=ll_chunksize),
                )
                .rename({"x": "lon", "y": "lat"})
                .squeeze("band", drop=True)
                .rename({"band_data": var})
                .drop_vars("spatial_ref")
                .expand_dims(**{key: [value] for key, value in meta.items()})
            )
            layer_ds.to_zarr(
                store_path, mode="a", region="auto", group=f"{dataset}/{var}"
            )


def setup_store(
    store_path,
    template_layer,
    dataset,
    var,
    grid_layers,
    metadata,
    ll_chunksize=8192,
):
    ds = xr.open_dataset(template_layer, engine="rasterio").rename(
        {"x": "lon", "y": "lat"}
    )
    logging.info("Opened template.")

    dims = ["lat", "lon"]
    coords = {
        "lat": list([float(l) for l in ds.coords["lat"]]),
        "lon": list([float(l) for l in ds.coords["lon"]]),
    }
    ds.close()
    del ds
    dim_lens = [len(coords["lat"]), len(coords["lon"])]
    chunk_sizes = [ll_chunksize, ll_chunksize]
    for dim in metadata["keys"]:
        dims.append(dim)
        coords[dim] = list(grid_layers[dim].unique())
        chunk_sizes.append(1)
        dim_lens.append(len(coords[dim]))

    logging.info("Got coords.")
    da = dask.array.zeros(
        dim_lens,
        chunks=chunk_sizes,
    )
    logging.info("Setup zero DataArray.")
    zero_ds = xr.Dataset(
        data_vars={var: (dims, da)},
        coords=coords,
    )
    logging.info("Composed zero dataset.")
    zero_ds.to_zarr(
        store_path,
        mode="w",
        group=f"{dataset}/{var}",
        safe_chunks=False,
        compute=False,
    )
    logging.info("Wrote zero dataset.")


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )
    main()
    logging.info("Done.")
