import logging
from fastapi import FastAPI
from fastapi.routing import APIRoute
from fastapi.logger import logger

from db import models
from db.database import engine

from .routers import attributes, features, tiles
from config import LOG_LEVEL

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


app = FastAPI(generate_unique_id_function=custom_generate_unique_id)

app.include_router(features.router, prefix="/features")
app.include_router(attributes.router, prefix="/attributes")
app.include_router(tiles.router, prefix="/tiles")
