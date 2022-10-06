"""
Global Configuration
"""

from os import getenv

RASTER_BASE_PATH = getenv("RASTER_BASE_PATH", "/data")
MYSQL_URI = getenv("MYSQL_URI")
API_TOKEN = getenv("API_TOKEN")
