"""Script to ingest generic rasters into terracotta.
"""

import csv
import json
from typing import Any, List
import argparse

import terracotta
import tqdm

parser = argparse.ArgumentParser(description="Terracotta Ingester")
parser.add_argument(
    "--csv",
    type=str,
    help="Path to the CSV file containing information for each raster",
)
parser.add_argument(
    "--metadata",
    type=str,
    help="Path to JSON file containing metadata",
)
parser.add_argument(
    "--local_path",
    type=str,
    help="Path to the raster file as far as this script is concerned",
)
parser.add_argument(
    "--db_path",
    type=str,
    help="Path to the raster file as the tileserver sees it, stored in terracotta database",
)


def read_csv(
    csv_filepath: str,
    local_path: str,
    db_path: str,
    metadata: dict,
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
                    "key_values": {_k: row[_k] for _k in metadata["variables"]},
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
                    f"Skipping {raster['path']} duplicate of existing dataset {dup_path}"
                )
                continue
            else:
                driver.insert(
                    raster["key_values"],
                    raster["local_path"],
                    override_path=raster["db_path"],
                )


if __name__ == "__main__":
    args = parser.parse_args()
    with open(args.metadata) as fh:
        metadata = json.load(fh)

    tile_keys: list[str] = metadata["keys"]
    raster_files = read_csv(
        args.input_csv_filepath,
        args.local_path,
        args.db_path,
    )
    ingest_files(metadata["source_db"], tile_keys, raster_files)
