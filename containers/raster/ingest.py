#!/usr/bin/env python3
import csv
import os
import sys
from typing import Dict, List
import argparse

import tqdm
import terracotta

parser = argparse.ArgumentParser(description="Terracotta Ingester")
parser.add_argument(
    "input_csv_filepath",
    type=str,
    help="Absolute path to the CSV file containing information for each raster",
)
parser.add_argument(
    "database_name",
    type=str,
    help="Name of the output database (will be created if it doesnt exist)",
)
parser.add_argument(
    "internal_raster_base_path", type=str, help="Internal path to the raster file"
)
parser.add_argument("keys", type=str, help="Internal path to the raster file")


# Define the location of the SQLite database
# (this will be created if it doesn't already exist)
# DB_NAME = f"{DB_PATH}/terracotta.sqlite"

# Define the list of keys that will be used to identify datasets.
# (these need to match the key_values dicts defined in RASTER_FILES below)
KEYS = ["type", "rp", "rcp", "epoch", "gcm"]


def load_csv(csv_filepath: str, internal_raster_base_path: str) -> List[dict]:
    """
    Define a list of raster files to import
    (this is a list of dictionaries, each with a file path and the values for
    each key - make sure the order matches the order of KEYS defined above)

    Read from hazard.csv
    hazard,path,rp,rcp,epoch,gcm,key
    """
    raster_files = []
    with open(csv_filepath) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            raster_files.append(
                {
                    "key_values": {
                        "type": row["hazard"],
                        "rp": row["rp"],
                        "rcp": row["rcp"],
                        "epoch": row["epoch"],
                        "gcm": row["gcm"],
                    },
                    "path": f"{internal_raster_base_path}/{row['key']}.tif",
                }
            )
    return raster_files


def ingest(db_name: str, keys: List[str], raster_files: List[Dict], append=True):
    driver = terracotta.get_driver(db_name)
    try:
        driver.create(keys)
    except Exception as err:
        print("Failed to create DB: %s", err)

    # sanity check that the database has the same keys that we want to load
    assert list(driver.key_names) == keys, (driver.key_names, keys)

    if append is True:
        # Remove existing rasters from the ingest list
        existing_rasters = driver.get_datasets().values()
        raster_files = [
            raster for raster in raster_files if raster["path"] not in existing_rasters
        ]

    progress_bar = tqdm.tqdm(raster_files)

    for idx, raster in enumerate(progress_bar):
        progress_bar.set_postfix(file=raster["path"])

        with driver.connect():
            try:
                category_map = {"meta_idx": str(idx)}
                metadata = driver.compute_metadata(
                    raster["path"], extra_metadata={"categories": category_map}
                )
                driver.insert(raster["key_values"], raster["path"], metadata=metadata)
            except Exception as err:
                print(
                    "raster {} failed ingest (skipping) due to {}".format(
                        raster["path"], err
                    )
                )


if __name__ == "__main__":
    args = parser.parse_args()
    raster_files = load_csv(args.input_csv_filepath, args.internal_raster_base_path)
    ingest(args.database_name, KEYS, raster_files)
