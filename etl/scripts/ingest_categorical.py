"""Script to ingest categorical rasters into terracotta.
"""

import csv
import json
import os
from typing import Any, List, OrderedDict
import argparse

import terracotta

parser = argparse.ArgumentParser(description="Terracotta Ingester")

parser.add_argument(
    "--raster",
    type=str,
    help="Path to a categorical raster file being loaded",
)
parser.add_argument(
    "--legend_csv",
    type=str,
    help="Absolute path to the CSV file containing legend information about a categorical raster being loaded",
)
parser.add_argument(
    "--metadata",
    type=str,
    help="Path to valid JSON file containing key:value mapping for the categorical raster.",
)
parser.add_argument(
    "--db_path",
    type=str,
    help="Path to the raster file as the tileserver sees it, stored in terracotta database",
)


def load_single_categorical(
    db_name: str,
    input_raster_fpath: str,
    db_path: str,
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
    tc_settings = terracotta.get_settings()
    tc_driver_path = f"{tc_settings.DRIVER_PATH}/{db_name}"
    driver = terracotta.get_driver(tc_driver_path, provider="postgresql")

    # Connect and setup DB
    driver.create(ordered_key_values.keys())

    # Load raster as per this: https://terracotta-python.readthedocs.io/en/latest/tutorials/categorical.html
    with driver.connect():
        metadata = driver.compute_metadata(
            input_raster_fpath, extra_metadata={"categories": category_map}
        )
    with driver.connect():
        db_filepath = f"{db_path}/{os.path.basename(input_raster_fpath)}"
        driver.insert(
            dict(ordered_key_values),
            input_raster_fpath,
            metadata=metadata,
            override_path=db_filepath,
        )


def read_legend_csv(fpath: str) -> dict:
    """
    Parse categorical values from a csv into a dictionary

    Assumes values are integers
    """
    _replaces = {"%": "perc", ">": "gt", "<": "lt"}
    output = {}
    with open(fpath, mode="r", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile)
        for line in reader:
            label = line["label"]
            for _find, _replace in _replaces.items():
                label = label.replace(_find, _replace)
            output[label] = int(line["value"])
    return output


if __name__ == "__main__":
    args = parser.parse_args()

    # Load and Parse the categorical CSV data
    categorical_map = read_legend_csv(
        args.legend_csv,
    )
    with open(args.metadata) as fh:
        metadata = json.load(fh)

    load_single_categorical(
        args.raster,
        args.db_path,
        metadata["keys"],
        categorical_map,
    )
