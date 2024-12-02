import logging

from fastapi import FastAPI
from fastapi.routing import APIRoute
from fastapi.logger import logger

import backend.app.routers.attributes as attributes
import backend.app.routers.features as features
import backend.app.routers.tiles as tiles
import backend.app.routers.colormap as colormap
from backend.config import LOG_LEVEL

formatter = logging.Formatter(
    "[%(asctime)s.%(msecs)03d] %(levelname)s %(filename)s - %(funcName)s - %(message)s",
    "%Y-%m-%d %H:%M:%S",
)
handler = logging.StreamHandler()
logging.getLogger().setLevel(LOG_LEVEL)
logger.addHandler(handler)
handler.setFormatter(formatter)


def custom_generate_unique_id(route: APIRoute):
    return f"{route.tags[0]}-{route.name}"


app = FastAPI(
    generate_unique_id_function=custom_generate_unique_id,
    title="GRI Infra-Risk-Vis API",
    description="API Supporting Global Resilience Initiative Visualisation UI.  Serving geospatial features (inc. related damages) and raster tiles (via Terracotta)",
    version="0.8.0",
    terms_of_service="",
    contact={
        "name": "Tom Russell",
        "url": "https://github.com/nismod/infra-risk-vis",
        "email": "",
    },
    license_info={
        "name": "MIT",
        "url": "https://raw.githubusercontent.com/nismod/infra-risk-vis/master/LICENSE",
    },
)

app.include_router(features.router, prefix="/features")
app.include_router(attributes.router, prefix="/attributes")
app.include_router(tiles.router, prefix="/tiles")
app.include_router(colormap.router, prefix="/colormap")
