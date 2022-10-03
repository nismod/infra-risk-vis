#!/usr/bin/env python3
import csv
import os
import sys
from typing import Dict, List

import tqdm
import terracotta

RASTER_BASE_PATH = sys.argv[1]
HAZARD_CSV = sys.argv[2]
DB_PATH = sys.argv[3]

# Define the location of the SQLite database
# (this will be created if it doesn't already exist)
DB_NAME = f"{DB_PATH}/terracotta.sqlite"

# Define the list of keys that will be used to identify datasets.
# (these need to match the key_values dicts defined in RASTER_FILES below)
KEYS = ["type", "rp", "rcp", "epoch", "gcm"]

# Define a list of raster files to import
# (this is a list of dictionaries, each with a file path and the values for
# each key - make sure the order matches the order of KEYS defined above)
#
# Read from hazard.csv
# hazard,path,rp,rcp,epoch,gcm,key
RASTER_FILES = []
with open(HAZARD_CSV) as fh:
    reader = csv.DictReader(fh)
    for row in reader:
        RASTER_FILES.append(
            {
                "key_values": {
                    "type": row["hazard"],
                    "rp": row["rp"],
                    "rcp": row["rcp"],
                    "epoch": row["epoch"],
                    "gcm": row["gcm"],
                },
                "path": f"{RASTER_BASE_PATH}/{row['key']}.tif",
            }
        )


def load(db_name: str, keys: List[str], raster_files: List[Dict], append=True):
    driver = terracotta.get_driver(db_name)

    # create an empty database if it doesn't exist
    if not os.path.isfile(db_name):
        driver.create(keys)

    # sanity check that the database has the same keys that we want to load
    assert list(driver.key_names) == keys, (driver.key_names, keys)

    if append is True:
        # Remove existing rasters from the ingest list
        existing_rasters = driver.get_datasets().values()
        raster_files = [
            raster for raster in raster_files if raster["path"] not in existing_rasters
        ]

    progress_bar = tqdm.tqdm(raster_files)

    for raster in progress_bar:
        progress_bar.set_postfix(file=raster["path"])

        with driver.connect():
            try:
                driver.insert(raster["key_values"], raster["path"])
            except Exception as err:
                print(
                    "raster {} failed ingest (skipping) due to {}".format(
                        raster["path"], err
                    )
                )


if __name__ == "__main__":
    load(DB_NAME, KEYS, RASTER_FILES)
