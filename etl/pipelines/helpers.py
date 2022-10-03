"""
Global Pipeline helpers
"""
import os
import logging


def get_logger(level=logging.DEBUG):
    # Create a custom logger
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(filename)s - %(funcName)s - %(levelname)s - %(message)s",
    )
    logger = logging.getLogger(__name__)
    return logger


def tiff_appears_valid(filepath: str) -> bool:
    """
    Check if the tiff file appears to be valid
    """
    if not os.path.exists(filepath):
        return False
    return os.path.getsize(filepath) > 0
