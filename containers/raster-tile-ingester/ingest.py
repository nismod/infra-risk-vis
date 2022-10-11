#!/usr/bin/env python3

"""
Script used for ingesting one or more rasters into a terracotta MySQL Database.

The run environment must contain:
    TC_DRIVER_PATH=mysql://USER:PASSWORD@HOST
"""
from copy import copy
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
    choices=["load_csv", "load_json", "load_single", "delete_database_entries"],
    help="Type of load operation - load_csv requires input_csv_filepath",
)
parser.add_argument(
    "--input_csv_filepath",
    type=str,
    help="Absolute path to the CSV file containing information for each raster",
)
parser.add_argument(
    "--tile_keys",
    type=str,
    help="A comma-seperated list of tile keys and ordering, e.g. hazard,rp,rcp,gcm",
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


def load_csv(
    csv_filepath: str,
    internal_raster_base_path: str,
    csv_key_column_map: dict,
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
        else:
            raise


def _build_mysql_driver_path(database: str, mysql_uri: str) -> str:
    """
    Build the full MySQL driver path for Terracotta using URI and database
    """
    return mysql_uri + "/" + database


def _build_sqlite_driver_path(database: str, path: str) -> str:
    """
    Build the full MySQL driver path for Terracotta using URI and database
    """
    return os.path.join(path, database + ".sqlite")


def _setup_driver(db_name: str) -> Any:
    """
    Set the terracotta driver
    """
    tc_settings = terracotta.get_settings()
    print(f"Using TC Settings: {tc_settings}")
    if tc_settings.DRIVER_PROVIDER == "mysql":
        tc_driver_path = _build_mysql_driver_path(db_name, tc_settings.DRIVER_PATH)
    else:
        tc_driver_path = _build_sqlite_driver_path(db_name, "/data/aqueduct")
    print(f"TC Driver Path: {tc_driver_path}")
    driver = terracotta.get_driver(tc_driver_path, provider=tc_settings.DRIVER_PROVIDER)
    return driver


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


def ingest_from_csv(
    db_name: str, keys: List[str], raster_files: List[dict], append=True
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

    progress_bar = tqdm.tqdm(raster_files)
    print("Starting ingest...")
    for raster in progress_bar:
        progress_bar.set_postfix(file=raster["path"])
        with driver.connect():
            try:
                # This does an internal DB check after each insert
                dup_path = _check_duplicate_entry(raster, driver)
                if dup_path:
                    print(
                        f"{raster['path']} has duplicate keys to existing dataset {dup_path} - will be skipped"
                    )
                    continue
                else:
                    driver.insert(raster["key_values"], raster["path"])
            except Exception as err:
                print(
                    "raster {} failed ingest (skipping) due to {}".format(
                        raster["path"], err
                    )
                )


def _delete_database_entries(db_name: str) -> bool:
    """
    Remove all datasets from the given database
    """
    print(f"Deleting all entries for {db_name}")
    driver = _setup_driver(db_name)
    count_before = len(driver.get_datasets())
    for ds_keys, ds_path in driver.get_datasets().items():
        driver.delete(ds_keys)
    count_after = len(driver.get_datasets())
    print(f"Deleted {count_before - count_after} entries from db {db_name}")


def _parse_csv_key_column_map(csv_key_column_map: str) -> dict:
    """
    Parse the key column map from args into a dict
    ::returns dict column_name_map e.g.:
        {
            "file_basename": "key",
            "type": "hazard",
            "rp": "rp",
            "rcp": "rcp",
            "epoch": "epoch",
            "gcm": "gcm",
        }
    """
    return json.loads(csv_key_column_map)


def _parse_tile_keys(tile_keys: str) -> List[str]:
    return tile_keys.split(",")


def _validate_keys_and_map(tile_keys: List[str], csv_key_column_map: dict):
    if not "file_basename" in csv_key_column_map.keys():
        raise Exception("csv_key_column_map must contain mapping for: 'file_basename'")
    _map = copy(csv_key_column_map)
    _map.pop("file_basename")
    if not sorted(tuple(tile_keys)) == sorted(_map.keys()):
        raise Exception("tile_keys do not match csv_key_column_map keys")


if __name__ == "__main__":
    args = parser.parse_args()
    if args.operation == "load_csv":
        csv_key_column_map = _parse_csv_key_column_map(args.csv_key_column_map)
        # Ensure keys in map match tile_keys
        tile_keys = _parse_tile_keys(args.tile_keys)
        _validate_keys_and_map(tile_keys, csv_key_column_map)
        print(f"parsed csv_key_column_map successfully as {csv_key_column_map}")
        raster_files = load_csv(
            args.input_csv_filepath, args.internal_raster_base_path, csv_key_column_map
        )
        ingest_from_csv(args.database_name, tile_keys, raster_files)
    elif args.operation == "delete_database_entries":
        _delete_database_entries(args.database_name)
    else:
        print(f"Operation {args.operation} not yet supported")
