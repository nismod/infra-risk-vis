"""
Global Pipeline helpers
"""
import csv
import os
import logging
from typing import List
from urllib.request import urlopen, urlretrieve


def get_logger(name: str = __name__, level=logging.INFO):
    # Create a custom logger
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(filename)s - %(funcName)s - %(levelname)s - %(message)s",
    )
    logger = logging.getLogger(name)
    return logger


def tiff_appears_valid(filepath: str) -> bool:
    """
    Check if the tiff file appears to be valid
    """
    if not os.path.exists(filepath):
        return False
    return os.path.getsize(filepath) > 0


def download_file(file_url: str, output_filepath: str) -> int:
    """
    Download the given file to given path

    ::return filesize int
    """
    urlretrieve(file_url, output_filepath)
    return os.path.getsize(output_filepath)


def file_key_from_fname(fname: str) -> str:
    """
    File name without file type
    """
    return fname.split(".")[0]


def fname_from_fpath(fpath: str) -> str:
    return os.path.basename(fpath)


def nc_fname_from_url(url: str) -> str:
    """
    Parse filename for nc file from url
    """
    return url.split("/")[-1]


def hazard_csv_exists(hazard_csv_fpath: str) -> bool:
    """
    Check the provided hazard csv exists
    """
    return os.path.exists(hazard_csv_fpath)


def hazard_csv_valid(hazard_csv_fpath: str, hazard_csv_fieldnames: List[str]) -> bool:
    """
    Ensure the provided hazard CSV is valid
    """
    with open(hazard_csv_fpath, "r") as csvfile:
        try:
            d_reader = csv.DictReader(csvfile)
            headers = d_reader.fieldnames
            assert set(headers) == set(
                hazard_csv_fieldnames
            ), "existing hazard csv has header mismatch: {} | {}".format(
                headers, hazard_csv_fieldnames
            )
        except:
            return False
        return True


def count_hazard_csv_rows(hazard_csv_fpath: str) -> int:
    """
    Count number of rows in the Hazard CSV
    """
    return sum(1 for _ in open(hazard_csv_fpath))
