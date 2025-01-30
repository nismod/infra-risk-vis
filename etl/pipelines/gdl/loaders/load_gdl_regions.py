"""
Load gdl_region DB table from JSON
"""

import sys
from typing import List
from loader_utils import load_json, init_db_session
from models import GdlRegion, IsoCountry


def load_gdl_regions(fpath: str) -> List:
    print("Loading with: ", fpath)
    data = load_json(fpath)
    engine, Session = init_db_session()

    loaded_ids = []
    with Session() as session:

        try:
            print("Deleting all rows in gdl_region table")
            session.query(GdlRegion).delete()

            iso_results = session.query(IsoCountry.iso_code)

            gdl_total_codes = []
            for result in iso_results:
                iso_code = result[0]
                gdl_total_code = iso_code + "t"
                gdl_total_codes.append(gdl_total_code)
                region = GdlRegion(
                    gdl_code=gdl_total_code,
                    region_name="Total",
                    level="National",
                    iso_code=iso_code,
                )
                session.add(region)
                _id = session.commit()
                loaded_ids.append(_id)

            for record in data:
                gdl_code = record["gdl_code"].lower()
                if gdl_code not in gdl_total_codes:
                    region = GdlRegion(
                        gdl_code=gdl_code,
                        region_name=record["region_name"],
                        level=record["level"],
                        iso_code=record["iso_code"].lower(),
                    )
                    session.add(region)
                    _id = session.commit()
                    loaded_ids.append(_id)

        except Exception as err:
            print(f"Region insert failed due to {err}, rolling back transaction...")
            session.rollback()
    engine.dispose()

    print(f"Loaded {len(loaded_ids)} regions to gdl_region")
    return loaded_ids


if __name__ == "__main__":
    if not len(sys.argv) == 2:
        print("Usage:", "load_gdl_regions.py <json file path>")
        sys.exit(1)

    fpath = sys.argv[1]
    if not fpath:
        print("missing fpath")
        sys.exit(1)

    loaded_ids = load_gdl_regions(fpath)
