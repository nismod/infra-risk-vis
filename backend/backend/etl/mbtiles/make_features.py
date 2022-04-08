from io import FileIO
import json
from typing import Iterable
from sqlalchemy.orm import Session
from geoalchemy2.shape import to_shape
from shapely.geometry import mapping

from backend.db.database import SessionLocal
from backend.db.models import Feature


def make_feature_properties(feature: Feature):
    properties = dict(feature.properties)
    for damage in feature.damages:
        key = f"{damage.hazard}__rcp_{damage.rcp}__epoch_{damage.epoch}__conf_None"
        properties[key] = damage.EAD_undefended_mean
    return properties


def feature_as_geojson(feature: Feature):
    properties = make_feature_properties(feature)
    return {
        "type": "Feature",
        "id": feature.id,
        "geometry": mapping(to_shape(feature.geom)),
        "properties": properties,
    }


def get_features_for_layer(layer: str):
    db: Session
    with SessionLocal() as db:
        for feature in db.query(Feature).filter(Feature.layer == layer):
            yield feature


def write_geojson_seq(geojsons: Iterable[dict], out: FileIO):
    for geojson in geojsons:
        out.write(json.dumps(geojson, indent=None) + "\n")
