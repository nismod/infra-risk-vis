#!/usr/bin/env python3
import csv
import json
import os
import sys
from typing import Any, Dict, List
import argparse
import traceback

import tqdm
import terracotta
from terracotta.exceptions import InvalidDatabaseError

parser = argparse.ArgumentParser(description="Terracotta Ingester")
parser.add_argument(
    "operation",
    type=str,
    choices=["load_csv", "load_json", "load_single"],
    help="Type of load operation - load_csv requires input_csv_filepath",
)
parser.add_argument(
    "--input_csv_filepath",
    type=str,
    help="Absolute path to the CSV file containing information for each raster",
)
parser.add_argument(
    "--csv_key_column_map",
    type=str,
    help='Map of DB Keys to column names for input CSV (must be valid JSON str and contain key: file_basename), e.g. \'{"file_basename": "key", "type": "hazard", "rp": "rp", "rcp": "rcp", "epoch": "epoch", "gcm": "gcm"}\'',
)
parser.add_argument(
    "--database_name",
    type=str,
    help="Name of the output database (will be created if it doesnt exist)",
)
parser.add_argument(
    "--internal_raster_base_path", type=str, help="Internal path to the raster file"
)


# Define the location of the SQLite database
# (this will be created if it doesn't already exist)
# DB_NAME = f"{DB_PATH}/terracotta.sqlite"

# Define the list of keys that will be used to identify datasets.
# (these need to match the key_values dicts defined in RASTER_FILES below)
# KEYS = ["type", "rp", "rcp", "epoch", "gcm"]


def build_driver_path(database: str, mysql_uri: str) -> str:
    """
    Build the full MySQL driver path for Terracotta using URI and database
    """
    return mysql_uri + "/" + database


def load_csv(
    csv_filepath: str,
    internal_raster_base_path: str,
    csv_key_column_map: dict = None,
) -> List[dict]:
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
                    # "key_values": {
                    #     "type": row["hazard"],
                    #     "rp": row["rp"],
                    #     "rcp": row["rcp"],
                    #     "epoch": row["epoch"],
                    #     "gcm": row["gcm"],
                    # },
                    "key_values": {
                        _k: row[_v]
                        for _k, _v in csv_key_column_map.items()
                        if _k != "file_basename"
                    },
                    "path": f"{internal_raster_base_path}/{row[csv_key_column_map['file_basename']]}.tif",
                }
            )
    return raster_files


def _create_db(db_name: str, driver: terracotta, keys: str) -> Any:
    """
    Create the Terracotta DB
    """
    try:
        print(f"Creating DB at path with keys {keys}")
        driver.create(keys)
        print(f"Create DB success")
    except InvalidDatabaseError as err:
        if "database exists" in traceback.format_exc():
            print(f"Database {db_name} already exists - skipped create")


def _setup_driver(db_name: str) -> Any:
    """
    Set the terracotta driver
    """
    tc_settings = terracotta.get_settings()
    tc_driver_path = build_driver_path(db_name, tc_settings.DRIVER_PATH)
    driver = terracotta.get_driver(tc_driver_path, provider=tc_settings.DRIVER_PROVIDER)
    return driver


def ingest_from_csv(
    db_name: str, keys: List[str], raster_files: List[Dict], append=True
):
    """
    Ingest to MySQL
    """
    # Setup Driver
    driver = _setup_driver(db_name)
    # Connect and setup DB
    _create_db(db_name, driver, keys)

    # sanity check that the database has the same keys that we want to load
    assert list(driver.key_names) == keys, (driver.key_names, keys)

    if append is True:
        # Remove existing rasters from the ingest list
        existing_rasters = driver.get_datasets().values()
        print(existing_rasters)
        raster_files = [
            raster for raster in raster_files if raster["path"] not in existing_rasters
        ]
        print(
            f"found {len(existing_rasters)} existing rasters in the DB - will append {len(raster_files)}"
        )

    progress_bar = tqdm.tqdm(raster_files[0:3])
    for idx, raster in enumerate(progress_bar):
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


def _parse_csv_key_colun_map(csv_key_colun_map: str) -> dict:
    """
    Parse the key column map from args into a dict
    """
    return json.loads(csv_key_colun_map)


if __name__ == "__main__":
    args = parser.parse_args()
    if args.operation == "load_csv":
        # csv_key_column_map = {
        #     "file_basename": "key",
        #     "type": "hazard",
        #     "rp": "rp",
        #     "rcp": "rcp",
        #     "epoch": "epoch",
        #     "gcm": "gcm",
        # }
        try:
            csv_key_column_map = _parse_csv_key_colun_map(args.csv_key_column_map)
            print(f"parsed csv_key_column_map successfully as {csv_key_column_map}")
        except:
            print(
                f"failed to parse csv_key_column_map, check if it is valid JSON: {args.csv_key_column_map}"
            )
        raster_files = load_csv(
            args.input_csv_filepath, args.internal_raster_base_path, csv_key_column_map
        )
        db_keys = list(csv_key_column_map.keys())
        db_keys.remove("file_basename")
        ingest_from_csv(args.database_name, db_keys, raster_files)
    else:
        print(f"Operation {args.operation} not yet supported")
