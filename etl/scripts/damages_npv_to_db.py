"""Load expected damages from a single file to database
"""
import pandas

from sqlalchemy.orm import Session
from tqdm import tqdm

from backend.db.database import SessionLocal
from backend.db.models import NPVDamage


def yield_npv_damages(data):
    """Generate NPVDamage database model objects."""
    for row in data.itertuples():
        yield NPVDamage(
            feature_id=row.uid,
            hazard=row.hazard,
            rcp=row.rcp,
            ead_amin=row.ead_amin,
            ead_mean=row.ead_mean,
            ead_amax=row.ead_amax,
            eael_amin=row.eael_amin,
            eael_mean=row.eael_mean,
            eael_amax=row.eael_amax,
        )


def parse_npv_damages(data):
    """Parse to tidy (long) dataframe"""
    # parquet files come in "fat" format, with string column names representing
    # hazard, RCP, damage type (EAD/EAEL), and min/mean/max suffix
    data_cols = [c for c in data.columns if "EA" in c]
    id_vars = ["uid"]

    # melt to long format
    data = (
        data.melt(id_vars=id_vars, value_vars=data_cols)
        .query("value > 0")
        .reset_index(drop=True)
    )

    # parse string key column for metadata
    meta = data.variable.str.extract(r"^([^_]+)__rcp_(\d.\d)__([^_]+)_(\w+)$")
    meta.columns = ["hazard", "rcp", "damage", "var"]

    # join metadata columns
    data = data.join(meta)

    # pivot back up so we end with a row per uid, hazard etc. (see index columns below)
    # and columns for each damage type, each with min/mean/max
    data = (
        data.drop(columns="variable")
        .pivot(
            index=["uid", "hazard", "rcp"],
            columns=["damage", "var"],
            values="value",
        )
        .fillna(0)
    )

    data.columns = [f"{var.lower()}_{stat}" for var, stat in data.columns]

    # ensure all columns are present - may be missing in case the data didn't
    # have any non-zero values in this batch
    expected_columns = [
        "ead_amin",
        "ead_mean",
        "ead_amax",
        "eael_amin",
        "eael_mean",
        "eael_amax",
    ]
    ensure_columns(data, expected_columns)

    return data.reset_index()


def ensure_columns(data, expected_columns):
    for col in expected_columns:
        if col not in data.columns:
            data[col] = 0
    return data


if __name__ == "__main__":
    try:
        input = snakemake.input
        npv_fname = snakemake.input.npv
        uid_fname = snakemake.input.uid
        layer_name = snakemake.wildcards.layer
        output = snakemake.output
        network_layers_fname = snakemake.config["network_layers"]
    except NameError:
        print("Expected to run from snakemake")
        exit()

    network_layers = pandas.read_csv(network_layers_fname)
    try:
        network_layer = network_layers[network_layers.damage_ref == layer_name].iloc[0]
    except IndexError as e:
        print(f"Could not find {layer_name} in network layers.")
        raise e

    damages = pandas.read_csv(npv_fname).set_index(network_layer.asset_id_column)
    ids = pandas.read_parquet(uid_fname).set_index(network_layer.asset_id_column)
    damages = damages.join(ids).reset_index(drop=True)

    damages_df = parse_npv_damages(damages)

    damages = yield_npv_damages(damages_df)

    db: Session
    with SessionLocal() as db:
        for i, damage_npv in enumerate(tqdm(damages, desc=f"{layer_name}_npv")):
            db.add(damage_npv)
            if i % 1000 == 0:
                db.commit()
        db.commit()

    with open(str(output), "w") as fh:
        fh.write(f"Loaded to database.\n\n")
        fh.write(f"From:\n{input}\n\n")
