"""
Global Configuration
"""

from os import getenv
import logging
import json

LOG_LEVEL = logging.getLevelName(getenv("LOG_LEVEL", "INFO"))
RASTER_BASE_PATH = getenv("RASTER_BASE_PATH", "/data")
MYSQL_URI = getenv("MYSQL_URI")
API_TOKEN = getenv("API_TOKEN")

# Temporarily here until we can re-eng the UI for domain / DB switching from API
# '{"fluvial": "aqueduct", "coastal": "aqueduct", "extreme_heat": "extreme_heat", "cyclone": "cyclone", "population": "jrc_pop"}'
DOMAIN_TO_DB_MAP = json.loads(getenv("DOMAIN_TO_DB_MAP", "{}"))
