"""
Load iso_country DB table from JSON
"""

import sys
from typing import List
from loader_utils import load_json, init_db_session
from models import IsoCountry


def load_iso_countries(fpath: str) -> List:
    print("Loading iso_country table from: ", fpath)
    data = load_json(fpath)
    engine, Session = init_db_session()

    loaded_ids = []
    with Session() as session:

        print("Deleting all rows in iso_country table")
        session.query(IsoCountry).delete()

        try:
            for record in data:
                country = IsoCountry(
                    iso_code=record["iso_code"].lower(),
                    country_name=record["country_name"],
                    continent=record["continent"],
                )
                session.add(country)
                _id = session.commit()
                loaded_ids.append(_id)
        except Exception as err:
            print(f"Country insert failed due to {err}, rolling back transaction...")
            session.rollback()
    engine.dispose()

    print(f"Loaded {len(loaded_ids)} countries to iso_country")
    return loaded_ids


if __name__ == "__main__":
    if not len(sys.argv) == 2:
        print("Usage:", "load_iso_countries.py <json file path>")
        sys.exit(1)

    fpath = sys.argv[1]
    if not fpath:
        print("missing fpath")
        sys.exit(1)

    loaded_ids = load_iso_countries(fpath)
