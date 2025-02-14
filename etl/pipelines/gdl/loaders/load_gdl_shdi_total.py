"""
Load GDL data from CSV
"""

from loader_utils import init_db_session
from models import IsoCountry, GdlRegion, GdlAnnual
import pandas as pd

import requests
from io import StringIO

DATA_URL = "https://zenodo.org/records/14868935/files/shdi_sgdi_total_8.0.csv"


def null_if_blank(value):
    return value if value != " " else None


def load_shdi_data():
    print("Loading GDL data from: ", DATA_URL)
    response = requests.get(DATA_URL)
    if response.status_code == 200:
        data = pd.read_csv(StringIO(response.text))
    else:
        print(f"Failed to fetch data. Status code: {response.status_code}")
        return

    loaded_country_ids = []
    loaded_region_ids = []
    loaded_record_ids = []

    engine, Session = init_db_session()
    with Session() as session:
        try:
            for _, row in data.iterrows():
                iso_code = row["iso_code"].lower()
                country_name = row["country"]
                continent = row["continent"]
                year = row["year"]
                gdl_code = row["GDLCODE"].lower()
                region_name = row["region"]
                level = row["level"]

                development = null_if_blank(row["shdi"])
                healthcare = null_if_blank(row["healthindex"])
                income = null_if_blank(row["incindex"])
                education = null_if_blank(row["edindex"])

                maybe_country = (
                    session.query(IsoCountry).filter_by(iso_code=iso_code).first()
                )
                if not maybe_country:
                    country = IsoCountry(
                        iso_code=iso_code,
                        country_name=country_name,
                        continent=continent,
                    )
                    session.add(country)
                    commit_id = session.commit()
                    loaded_country_ids.append(commit_id)

                maybe_region = (
                    session.query(GdlRegion).filter_by(gdl_code=gdl_code).first()
                )
                if not maybe_region:
                    region = GdlRegion(
                        gdl_code=gdl_code,
                        region_name=region_name,
                        level=level,
                        iso_code=iso_code,
                    )
                    session.add(region)
                    commit_id = session.commit()
                    loaded_region_ids.append(commit_id)

                year_record = GdlAnnual(
                    gdl_code=gdl_code,
                    year=year,
                    development=development,
                    education=education,
                    healthcare=healthcare,
                    income=income,
                )
                session.add(year_record)
                commit_id = session.commit()
                loaded_record_ids.append(commit_id)

        except Exception as err:
            print(f"Data insert failed due to {err}, rolling back transaction...")
            session.rollback()

        engine.dispose()
        print(f"Loaded {len(loaded_country_ids)} ISO country rows")
        print(f"Loaded {len(loaded_region_ids)} GDL region rows")
        print(f"Loaded {len(loaded_record_ids)} annual GDL data rows")


if __name__ == "__main__":
    load_shdi_data()
