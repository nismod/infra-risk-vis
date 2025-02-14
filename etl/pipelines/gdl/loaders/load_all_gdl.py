"""
Load GDL annual data, national boundaries, and subnational boundaries
"""

from load_gdl_national import load_gdl_national
from load_gdl_subnational import load_gdl_subnational
from loader_utils import init_db_session
from models import (
    GdlAnnual,
    GdlNational,
    GdlSubnational,
    GdlRegion,
    IsoCountry,
)

from load_gdl_shdi_total import load_shdi_data


def load_all_data():
    engine, Session = init_db_session()

    # Delete everything to start
    with Session() as session:
        session.query(GdlAnnual).delete()
        session.query(GdlNational).delete()
        session.query(GdlSubnational).delete()
        session.query(GdlRegion).delete()
        session.query(IsoCountry).delete()
        session.commit()

    engine.dispose()

    load_shdi_data()
    load_gdl_subnational()
    load_gdl_national()


if __name__ == "__main__":
    load_all_data()
