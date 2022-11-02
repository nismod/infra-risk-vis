"""Write from database to GeoJSONSeq file for a single visualisation layer
"""
import os
import sys
from operator import attrgetter

import ujson as json

from geoalchemy2.shape import to_shape
from shapely.geometry import mapping
from sqlalchemy.orm import Session, selectinload
from tqdm import tqdm

from backend.db.database import SessionLocal
from backend.db.models import Feature


def feature_as_geojson(feature: Feature):
    properties = {
        "asset_id": feature.properties["asset_id"],
        "asset_type": feature.properties["asset_type"],
    }
    # sort so any higher protection standard will overwrite
    damages = sorted(feature.damages_expected, key=attrgetter("protection_standard"))
    for damage in damages:
        key = f"{damage.hazard}__rcp_{damage.rcp}__epoch_{damage.epoch}__conf_None"
        properties[f"ead__{key}"] = damage.ead_mean
        properties[f"eael__{key}"] = damage.eael_mean

    return {
        "type": "Feature",
        "id": feature.id,
        "geometry": mapping(to_shape(feature.geom)),
        "properties": properties,
    }


def yield_features_for_layer(layer: str):
    db: Session
    with SessionLocal() as db:
        query = (
            db.query(Feature)
            .options(selectinload(Feature.damages_expected))
            .filter(Feature.layer == layer)
            .execution_options(yield_per=1000)
        )
        for partition in db.execute(query).partitions(1000):
            for (feature,) in partition:
                yield feature


if __name__ == "__main__":
    try:
        layer = snakemake.wildcards.layer
        output = snakemake.output
    except NameError:
        print("Expected to run from snakemake")
        exit()

    with open(str(output), "w") as fh:
        for feature in tqdm(yield_features_for_layer(layer)):
            geojson = feature_as_geojson(feature)
            json.dump(geojson, fh, indent=0, ensure_ascii=False)
            fh.write("\n")
