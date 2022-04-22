"""Load feature layer definitions to database.
"""
import pandas

from sqlalchemy.orm import Session

from backend.db.database import SessionLocal
from backend.db.models import FeatureLayer


if __name__ == "__main__":
    try:
        output = snakemake.output
        network_tilelayers_fname = snakemake.config["network_tilelayers"]
    except NameError:
        print("Expected to run from snakemake")
        exit()

    network_tilelayers = pandas.read_csv(network_tilelayers_fname)
    db: Session
    with SessionLocal() as db:
        for row in network_tilelayers.itertuples():
            feature_layer = FeatureLayer(
                layer_name=row.layer,
                sector=row.sector,
                subsector=row.subsector,
                asset_type=row.asset_type,
            )
            db.add(feature_layer)
        db.commit()

    with open(str(output), "w") as fh:
        fh.write(f"Loaded to database.\n\n")
        fh.write(
            f"From:\n{network_tilelayers_fname}\n\n"
        )
