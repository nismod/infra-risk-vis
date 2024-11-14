"""
Singleton Helpers
"""

import logging
import traceback

from config import TILEDB_URI


def build_driver_path(database: str, tiledb_uri: str = TILEDB_URI) -> str:
    """
    Build the full MySQL driver path for Terracotta using URI and database
    """
    return tiledb_uri + "/" + database


def handle_exception(logger: logging.Logger, err: Exception):
    """
    Handle generic exceptions
    """
    logger.error(
        "%s failed with tb %s, error: %s",
        __name__,
        traceback.format_tb(err.__traceback__),
        err,
    )
