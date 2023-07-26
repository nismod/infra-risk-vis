#!/usr/bin/env python3

"""
Script used for ingesting one or more rasters into a terracotta MySQL Database.
"""
from copy import copy
import csv
import json
import os
import sys
from typing import Any, Dict, List, OrderedDict, Sequence
import argparse
import traceback
import pymysql

import tqdm
import terracotta
from terracotta.exceptions import InvalidDatabaseError

parser = argparse.ArgumentParser(description="Terracotta Ingester")
parser.add_argument(
    "operation",
    type=str,
    choices=[
        "load_csv",
        "load_single_categorical",
        "delete_database_entries",
        "drop_database",
    ],
    help="""Type of load operation - 
        - `load_csv` Loads all rasters in a CSV file generated from an ETL pipeline, e.g. with the header:  `hazard,metric,path,rcp,epoch,gcm,key`
        - `load_single_categorical` Loads a single categorical raster (using provided category_map).  Requires: `input_raster_filepath`, `categorical_legend_csv_filepath`, `categorical_csv_label_column`, `categorical_csv_value_column`
        - `delete_database_entries` Delete all raster entries from a given database (leaving the database empty).  Requires `database_name`
        - `drop_database` Drop a database and all raster entries contained-within.  Requires `database_name`
    """,
)
parser.add_argument(
    "--input_csv_filepath",
    type=str,
    help="Absolute path to the CSV file containing information for each raster",
)
parser.add_argument(
    "--input_raster_filepath",
    type=str,
    help="Absolute path to a categorical raster file being loaded",
)
parser.add_argument(
    "--categorical_legend_csv_filepath",
    type=str,
    help="Absolute path to the CSV file containing legend information about a categorical raster being loaded",
)
parser.add_argument(
    "--categorical_csv_label_column",
    type=str,
    help="Name of the label column in CSV",
)
parser.add_argument(
    "--categorical_csv_value_column",
    type=str,
    help="Name of the value column in CSV",
)
parser.add_argument(
    "--categorical_key_values_json_path",
    type=str,
    help='Path to valid JSON file containing key:value mapping for the categorical raster.  NOTE: this will be loaded as an OrderedDict - so the key ordering will be maintained for DB usage.',
)
parser.add_argument(
    "--tile_keys_path",
    type=str,
    help='Path to JSON file containing list of tile keys and ordering, e.g. ["hazard", "rp", "rcp", gcm"]',
)
parser.add_argument(
    "--csv_to_db_field_map_path",
    type=str,
    help=(
        'Path to JSON file with map of DB Keys to metadata.CSV column names (must contain keys: "file_basename" and "type"), e.g. '
        '{"file_basename": "key", "type": "hazard", "rp": "rp", "rcp": "rcp", "epoch": "epoch", "gcm": "gcm"} '
    )
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
                    "key_values": {
                        _k: row[_v]
                        for _k, _v in csv_key_column_map.items()
                        if _k != "file_basename"
                    },
                    "path": f"{internal_raster_base_path}/{row[csv_key_column_map['file_basename']]}.tif",
                }
            )
    return raster_files


def _drop_db(db_name: str) -> None:
    """
    Drop the given database - must have root perms
    """
    import pymysql

    tc_settings = terracotta.get_settings()
    mysql_uri = tc_settings.DRIVER_PATH
    if not "mysql://" in mysql_uri:
        raise Exception("can only drop databaess from mysql")
    parts = mysql_uri.replace("mysql://", "").split("@")
    host = parts[1]
    username, password = parts[0].split(":")
    # Connect to the database
    connection = pymysql.connect(
        host=host,
        user=username,
        password=password,
        database="mysql",
        charset="utf8mb4",
        cursorclass=pymysql.cursors.DictCursor,
    )

    try:
        with connection:
            with connection.cursor() as cursor:
                # Create a new record
                sql = f"DROP DATABASE {db_name}"
                cursor.execute(sql)
            connection.commit()
            print(f"{db_name=} deleted.")
    except pymysql.err.OperationalError as e:
        if str(pymysql.constants.ER.DB_DROP_EXISTS) in repr(e):
            print(f"{db_name=} did not exist, couldn\'t delete.")
        else:
            raise e


def _create_db(db_name: str, driver: terracotta, keys: Sequence[str]) -> Any:
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


def _check_driver_keys(driver: Any, keys: List[str]) -> None:
    """
    Check the driver keys match the given keys
    """
    assert sorted(list(driver.key_names)) == sorted(keys), (driver.key_names, keys)


def _load_single_categorical(
    db_name: str,
    input_raster_fpath: str,
    internal_raster_base_path: str,
    ordered_key_values: OrderedDict,
    category_map: dict,
) -> None:
    """
    Load a single categorical raster

    ::param db_name str Database name
    ::param input_raster_fpath str Filepath to input categorical raster
    ::param internal_raster_base_path str The Internal path the raster will be stored-under
    ::param ordered_key_values: OrderedDict e.g. {'type' : 'exposure', 'sensor':'S2', 'data':'20181010', 'band':'cloudmask'}
    ::param category_map dict The category map (label:category_value) for the raster e.g.:
        {
            'clear land': 0,
            'clear water': 1,
            'cloud': 2,
            'cloud shadow': 3
        }
    """
    # Setup Driver
    driver = _setup_driver(db_name)
    # Connect and setup DB
    _create_db(db_name, driver, ordered_key_values.keys())
    # sanity check that the database has the same keys that we want to load
    _check_driver_keys(driver, ordered_key_values.keys())
    # Load raster as per this: https://terracotta-python.readthedocs.io/en/latest/tutorials/categorical.html
    with driver.connect():
        print(
            f"Computing metadata for categorical raster: {os.path.basename(input_raster_fpath)}"
        )
        metadata = driver.compute_metadata(
            input_raster_fpath, extra_metadata={"categories": category_map}
        )
    with driver.connect():
        internal_raster_fpath = (
            f"{internal_raster_base_path}/{os.path.basename(input_raster_fpath)}",
        )
        print(f"Inserting categorical raster: {os.path.basename(input_raster_fpath)}")
        driver.insert(
            dict(ordered_key_values), internal_raster_fpath, metadata=metadata
        )
        print(f"Inserted: {os.path.basename(input_raster_fpath)}")
    print("Completed categorical single ingest")


def _parse_categorical_csv(
    cat_csv_fpath: str,
    label_column: str,
    value_column: str,
) -> dict:
    """
    Parse categorical values from a csv into a dictionary

    Assumes values are integers
    """
    _replaces = {"%": "perc", ">": "gt", "<": "lt"}
    output = {}
    with open(cat_csv_fpath, mode="r", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile)
        for line in reader:
            label = line[label_column]
            for _find, _replace in _replaces.items():
                label = label.replace(_find, _replace)
            output[label] = int(line[value_column])
    return output


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
    _check_driver_keys(driver, keys)

    progress_bar = tqdm.tqdm(raster_files)
    print("Starting ingest...")
    with driver.connect():
        for raster in progress_bar:
            progress_bar.set_postfix(file=raster["path"])
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


def _parse_csv_key_column_map(csv_key_column_map_path: str) -> dict:
    """
    Parse the key column map file from args into a dict
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
    with open(csv_key_column_map_path) as fp:
        data = json.load(fp)
    return data


def _parse_tile_keys(tile_key_path: str) -> List[str]:
    with open(tile_key_path) as fp:
        data = json.load(fp)
    return data


def _parse_ordered_key_values(json_filepath: str) -> OrderedDict:
    """
    Parse JSON file into Ordered Dict
    """
    import collections

    with open(json_filepath, "r") as fp:
        data = json.load(fp, object_pairs_hook=collections.OrderedDict)

    return data


def _validate_keys_and_map(tile_keys: List[str], csv_key_column_map: dict):
    if not "file_basename" in csv_key_column_map.keys():
        raise Exception("DB field to CSV column name mapping requires 'file_basename' key")
    if not "type" in csv_key_column_map.keys():
        raise Exception("DB field to CSV column name mapping requires 'type' key")
    _map = copy(csv_key_column_map)
    _map.pop("file_basename")
    if not sorted(tuple(tile_keys)) == sorted(_map.keys()):
        raise Exception("tile keys do not match DB field to CSV column keys")


if __name__ == "__main__":
    args = parser.parse_args()
    if args.operation == "load_csv":
        csv_key_column_map: dict[str, str] = _parse_csv_key_column_map(args.csv_to_db_field_map_path)
        # Ensure keys in map match tile_keys
        tile_keys: list[str] = _parse_tile_keys(args.tile_keys_path)
        _validate_keys_and_map(tile_keys, csv_key_column_map)
        print(f"parsed csv_key_column_map successfully as {csv_key_column_map}")
        raster_files = load_csv(
            args.input_csv_filepath, args.internal_raster_base_path, csv_key_column_map
        )
        ingest_from_csv(args.database_name, tile_keys, raster_files)
    elif args.operation == "delete_database_entries":
        _delete_database_entries(args.database_name)
    elif args.operation == "drop_database":
        _drop_db(args.database_name)
    elif args.operation == "load_single_categorical":
        args_required = [
            args.categorical_legend_csv_filepath,
            args.categorical_csv_label_column,
            args.categorical_csv_value_column,
            args.database_name,
            args.input_raster_filepath,
            args.internal_raster_base_path,
            args.categorical_key_values_json_path,
        ]
        if None in args_required:
            print(
                """Following args must all be set for categorical load: 
                --categorical_legend_csv_filepath,
                --categorical_csv_label_column,
                --categorical_csv_value_column,
                --database_name,
                --input_raster_filepath,
                --internal_raster_base_path,
                --categorical_key_values_json_path"""
            )
        else:
            # Load and Parse the categorical CSV data
            categorical_map = _parse_categorical_csv(
                args.categorical_legend_csv_filepath,
                args.categorical_csv_label_column,
                args.categorical_csv_value_column,
            )
            ordered_key_values = _parse_ordered_key_values(args.categorical_key_values_json_path)
            _load_single_categorical(
                args.database_name,
                args.input_raster_filepath,
                args.internal_raster_base_path,
                ordered_key_values,
                categorical_map,
            )
    else:
        print(f"Operation {args.operation} not yet supported")
    print("Done.")
