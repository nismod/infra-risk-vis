"""
Load gdl_national DB table from GeoJSON
"""

import json
from typing import List
import sqlalchemy as sa
from loader_utils import load_json, init_db_session
from models import GdlNational, GdlRegion

import requests

DATA_URL = "https://zenodo.org/records/14868935/files/gdl_v6.4_national_small.json"


def load_gdl_national() -> List:
    print("Loading GDL data from: ", DATA_URL)
    response = requests.get(DATA_URL)
    if response.status_code == 200:
        data = response.json()
    else:
        print(f"Failed to fetch data. Status code: {response.status_code}")
        return

    loaded_ids = []
    engine, Session = init_db_session()
    with Session() as session:
        try:
            for feature in data["features"]:
                synthetic_gdl_code = (feature["properties"]["iso_code"] + "t").lower()

                maybe_region = (
                    session.query(GdlRegion)
                    .filter_by(gdl_code=synthetic_gdl_code)
                    .first()
                )

                if maybe_region:
                    boundary = GdlNational(
                        gdl_code=synthetic_gdl_code,
                        # format geojson for PostGIS
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
                    commit_id = session.commit()
                    loaded_ids.append(commit_id)

                else:
                    print(
                        f"Skipped (due to no matching gdl_code in GdlRegion): gdl_code:'{synthetic_gdl_code}'"
                    )

        except Exception as err:
            print(f"Boundary insert failed due to {err}, rolling back transaction...")
            session.rollback()
    engine.dispose()

    print(f"Loaded {len(loaded_ids)} boundaries to gdl_national")


if __name__ == "__main__":
    load_gdl_national()
