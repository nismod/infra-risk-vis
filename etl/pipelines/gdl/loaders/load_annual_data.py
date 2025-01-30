"""
Load GDL annual data from JSONs
"""

import sys
from typing import List
from loader_utils import load_json, init_db_session
from models import GdlAnnual, GdlRegion


def add_years_to_set(data, included_years: set):
    for record in data:
        for key in record:
            if key.isnumeric():
                included_years.add(key)


def add_region_years_to_set(data, included_years: set, gdl_code: str):
    for record in data:
        if record["GDLCODE"].lower() == gdl_code:
            for key in record:
                if key.isnumeric():
                    included_years.add(key)


def add_gdl_to_set(data, included_gdl: set):
    for record in data:
        gdl_code = record["GDLCODE"].lower()
        included_gdl.add(gdl_code)


def load_dataset(session, Model, dataset, dataset_key, loaded_ids):
    for record in dataset:
        gdl_code = record["GDLCODE"].lower()
        region_result = (
            session.query(GdlRegion).where(GdlRegion.gdl_code == gdl_code).first()
        )
        if not region_result:
            print("No matching gdl_code:", gdl_code)

        else:
            for key in record:
                # GDL download is denormalized (entry per GDL region code), numeric keys IFF they pertain to annual vals
                if key.isnumeric():
                    value = record[key]
                    if value is not None:
                        annual_result = (
                            session.query(Model)
                            .where(Model.gdl_code == gdl_code)
                            .where(Model.year == key)
                            .first()
                        )

                        if annual_result is not None:
                            # update existing row
                            session.query(Model).where(
                                Model.gdl_code == gdl_code
                            ).where(Model.year == key).update({dataset_key: value})

                        else:
                            # add new row
                            entry = Model(
                                gdl_code=record["GDLCODE"].lower(),
                                year=key,
                                **{dataset_key: value},
                            )
                            session.add(entry)

                        _id = session.commit()
                        loaded_ids.append(_id)


def load_annual_data(
    development_fpath: str,
    education_fpath: str,
    income_fpath: str,
    healthcare_fpath: str,
    Model,
) -> List:
    development_data = load_json(development_fpath)
    education_data = load_json(education_fpath)
    income_data = load_json(income_fpath)
    healthcare_data = load_json(healthcare_fpath)

    datasets = [
        {"data": development_data, "key": "development"},
        {"data": education_data, "key": "education"},
        {"data": income_data, "key": "income"},
        {"data": healthcare_data, "key": "healthcare"},
    ]

    engine, Session = init_db_session()

    loaded_ids = []
    with Session() as session:
        try:
            print("Deleting all rows in table:", Model.__table__.name)
            session.query(Model).delete()

            for dataset in datasets:
                load_dataset(
                    session=session,
                    Model=Model,
                    dataset=dataset["data"],
                    dataset_key=dataset["key"],
                    loaded_ids=loaded_ids,
                )

        except Exception as err:
            print(f"Data entry insert failed due to {err}, rolling back transaction...")
            session.rollback()

    engine.dispose()

    print(f"Loaded {len(loaded_ids)} records to table:", Model.__table__.name)
    return loaded_ids


def load_all_tables(development_fpath, education_fpath, income_fpath, healthcare_fpath):
    load_annual_data(
        development_fpath, education_fpath, income_fpath, healthcare_fpath, GdlAnnual
    )


if __name__ == "__main__":
    if not len(sys.argv) == 5:
        print(
            "Usage:",
            "load_annual_data.py <development json> <education json> <income json> <healthcare json>",
        )
        sys.exit(1)

    fpath = sys.argv[1]
    if not fpath:
        print("missing fpath")
        sys.exit(1)

    load_all_tables(
        development_fpath=sys.argv[2],
        education_fpath=sys.argv[3],
        income_fpath=sys.argv[4],
        healthcare_fpath=sys.argv[5],
    )
