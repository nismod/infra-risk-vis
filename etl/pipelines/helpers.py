"""
Global Pipeline helpers
"""
import os
import logging
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


def _download_file(file_url: str, output_filepath: str) -> int:
    """
    Download the given file to given path

    ::return filesize int
    """
    urlretrieve(file_url, output_filepath)
    return os.path.getsize(output_filepath)
