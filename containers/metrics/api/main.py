import logging

from fastapi import FastAPI
from fastapi.routing import APIRoute
from fastapi.logger import logger
from fastapi.middleware.cors import CORSMiddleware

from api.routers.gdl import router

formatter = logging.Formatter(
    "[%(asctime)s.%(msecs)03d] %(levelname)s %(filename)s - %(funcName)s - %(message)s",
    "%Y-%m-%d %H:%M:%S",
)
handler = logging.StreamHandler()
logging.getLogger().setLevel("DEBUG")
logger.addHandler(handler)
handler.setFormatter(formatter)


def custom_generate_unique_id(route: APIRoute):
    return f"{route.tags[0]}-{route.name}"


app = FastAPI(
    generate_unique_id_function=custom_generate_unique_id,
    title="GRI Metrics API",
    description="API Supporting Global Resilience Initiative country metrics UI.  Serving geojson features and json metrics",
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

origins = [
    "http://localhost",
    "http://localhost:5173",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(router)
