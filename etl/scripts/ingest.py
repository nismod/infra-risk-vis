"""Script to ingest generic rasters into terracotta.
"""

import csv
import json
from typing import Any, List

import terracotta
import tqdm


def read_csv(
    csv_filepath: str,
    tile_keys: list[str],
    local_path: str,
    db_path: str,
) -> List[dict]:
    """
    Define a list of raster files to import

    Read from hazard.csv
    hazard,path,rp,rcp,epoch,gcm,filename
    """
    raster_files = []
    with open(csv_filepath) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            raster_files.append(
                {
                    "key_values": {_k: row[_k] for _k in tile_keys},
                    "local_path": f"{local_path}/{row['filename']}",
                    "db_path": f"{db_path}/{row['filename']}",
                }
            )
    return raster_files


def _check_duplicate_entry(raster_file: dict, driver: Any) -> str:
    """
    Check if the given raster file has dup keys to something already in the DB

    terracotta existing rasters:
        e.g.: {('coastal', '1', 'baseline', 'present', 'None'): '/data/aqueduct/inuncoast_historical_nosub_hist_rp0001_5.tif'
    """
    for ds_keys, ds_path in driver.get_datasets().items():
        if sorted(tuple(raster_file["key_values"].values())) == sorted(ds_keys):
            return ds_path
    return None


def ingest_files(db_name: str, keys: List[str], raster_files: List[dict]):
    # Setup Driver
    tc_settings = terracotta.get_settings()
    tc_driver_path = f"{tc_settings.DRIVER_PATH}/{db_name}"
    driver = terracotta.get_driver(tc_driver_path, "postgresql")

    # Connect and setup DB
    driver.create(keys)

    progress_bar = tqdm.tqdm(raster_files)
    with driver.connect():
        for raster in progress_bar:
            progress_bar.set_postfix(file=raster["local_path"])

            # This does an internal DB check after each insert
            dup_path = _check_duplicate_entry(raster, driver)
            if dup_path:
                print(
                    f"Skipping {raster['local_path']} duplicate of existing dataset {dup_path}"
                )
                continue
            else:
                driver.insert(
                    raster["key_values"],
                    raster["local_path"],
                    override_path=raster["db_path"],
                )


if __name__ == "__main__":
    try:
        metadata_path = snakemake.input.metadata
        layers_path = snakemake.input.layers
        flag = snakemake.output.flag
    except NameError:
        assert False, "Must be run from snakemake"

    with open(metadata_path) as fh:
        metadata = json.load(fh)

    tile_keys: list[str] = metadata["keys"]
    raster_files = read_csv(
        layers_path,
        tile_keys,
        f"raster/cog/{metadata['domain']}",
        f"/data/{metadata['domain']}",
    )
    ingest_files(metadata["source_db"], tile_keys, raster_files)

    with open(flag, "w") as fh:
        fh.write("Done")
