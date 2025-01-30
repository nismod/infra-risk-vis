import traceback

import os
from dotenv import load_dotenv

load_dotenv()


def get_log_level():
    """
    Return the .env log level or default
    """
    maybe_log_level = os.getenv("METRICS_LOG_LEVEL")
    return maybe_log_level if maybe_log_level != None else "DEBUG"


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
