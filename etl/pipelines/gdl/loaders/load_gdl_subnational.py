"""
Load gdl_subnational DB table from GeoJSON
"""

import sys
import json
from typing import List
import sqlalchemy as sa
from loader_utils import load_json, init_db_session
from models import GdlSubnational


def load_gdl_subnational(fpath: str) -> List:
    print("Loading gdl_subnational from: ", fpath)
    data = load_json(fpath)
    engine, Session = init_db_session()

    loaded_ids = []
    with Session() as session:

        print("Deleting all rows in gdl_subnational table")
        session.query(GdlSubnational).delete()

        try:
            for feature in data["features"]:
                boundary = GdlSubnational(
                    gdl_code=feature["properties"]["gdlcode"].lower(),
                    geometry=sa.func.ST_AsEWKT(
                        sa.func.ST_SetSRID(
                            sa.func.ST_Multi(
                                sa.func.ST_GeomFromGeoJSON(
                                    json.dumps(feature["geometry"])
                                )
                            ),
                            4326,
                        )
                    ),
                )
                session.add(boundary)
                _id = session.commit()
                loaded_ids.append(_id)
        except Exception as err:
            print(f"Boundary insert failed due to {err}, rolling back transaction...")
            session.rollback()
    engine.dispose()

    print(f"Loaded {len(loaded_ids)} boundaries to gdl_subnational")
    return loaded_ids


if __name__ == "__main__":
    if not len(sys.argv) == 2:
        print("Usage:", "load_gdl_subnational.py <geojson file path>")
        sys.exit(1)

    fpath = sys.argv[1]
    if not fpath:
        print("missing fpath")
        sys.exit(1)

    loaded_ids = load_gdl_subnational(fpath)
