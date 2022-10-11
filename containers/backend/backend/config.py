"""
Global Configuration
"""

from os import getenv

RASTER_BASE_PATH = getenv("RASTER_BASE_PATH", "/data")
MYSQL_URI = getenv("MYSQL_URI")
API_TOKEN = getenv("API_TOKEN")

# Temporarily here until we can re-eng the UI for domain / DB switching from API
DOMAIN_TO_DB_MAP = {
    "fluvial": "aqueduct",
    "coastal": "aqueduct",
    "extreme_heat": "extreme_heat",
    "cyclone": "cyclone",
}
