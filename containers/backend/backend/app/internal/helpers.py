"""
Singleton Helpers
"""
import traceback

from config import MYSQL_URI


def build_driver_path(database: str, mysql_uri: str = MYSQL_URI) -> str:
    """
    Build the full MySQL driver path for Terracotta using URI and database
    """
    return mysql_uri + "/" + database


def handle_exception(logger, err: Exception):
    """
    Handle generic exceptions
    """
    logger.error(
        "%s failed with tb %s, error: %s",
        __name__,
        traceback.format_tb(err.__traceback__),
        err,
    )
