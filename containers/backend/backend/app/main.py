import logging

from fastapi import FastAPI
from fastapi.routing import APIRoute
from fastapi.logger import logger

from db import models
from db.database import engine
from .routers import attributes, features, tiles, colormap
from config import LOG_LEVEL, DOMAIN_TO_DB_MAP

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

@app.on_event("startup")
async def log_domain_to_db_mapping():
    logging.info(f"Using frontend 'domain' to tiles-db database name mapping:\n{DOMAIN_TO_DB_MAP}")

app.include_router(features.router, prefix="/features")
app.include_router(attributes.router, prefix="/attributes")
app.include_router(tiles.router, prefix="/tiles")
app.include_router(colormap.router, prefix="/colormap")
