"""
Load GDL annual data from JSON
"""

from load_iso_countries import load_iso_countries
from load_gdl_regions import load_gdl_regions
from load_gdl_national import load_gdl_national
from load_gdl_subnational import load_gdl_subnational
from load_annual_data import load_all_tables
from loader_utils import init_db_session
from models import (
    GdlAnnual,
    GdlNational,
    GdlSubnational,
    GdlRegion,
    IsoCountry,
)

# Data paths
ISO_META_PATH = "../json/iso_countries.json"
REGIONS_META_PATH = "../json/gdl_regions.json"
NATIONAL_GEOJSON_PATH = "../geojson/GDL_V6.3_national_small.json"
SUBNATIONAL_GEOJSON_PATH = "../geojson/gdl_v6.3_large_visvaligram_weighted_0.02.json"
DEVELOPMENT_ANNUAL_PATH = "../json/development.json"
EDUCATION_ANNUAL_PATH = "../json/education.json"
HEALTHCARE_ANNUAL_PATH = "../json/healthcare.json"
INCOME_ANNUAL_PATH = "../json/income.json"


def load_all_data():
    engine, Session = init_db_session()

    with Session() as session:
        session.query(GdlAnnual).delete()
        session.query(GdlNational).delete()
        session.query(GdlSubnational).delete()
        session.query(GdlRegion).delete()
        session.query(IsoCountry).delete()
        session.commit()

    engine.dispose()

    load_iso_countries(ISO_META_PATH)
    load_gdl_regions(REGIONS_META_PATH)
    load_gdl_subnational(SUBNATIONAL_GEOJSON_PATH)
    load_gdl_national(NATIONAL_GEOJSON_PATH)
    load_all_tables(
        development_fpath=DEVELOPMENT_ANNUAL_PATH,
        education_fpath=EDUCATION_ANNUAL_PATH,
        income_fpath=INCOME_ANNUAL_PATH,
        healthcare_fpath=HEALTHCARE_ANNUAL_PATH,
    )


if __name__ == "__main__":
    load_all_data()
