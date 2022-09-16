"""Load expected damages from a single file to database
"""
import pandas

from sqlalchemy.orm import Session
from tqdm import tqdm

from backend.db.database import SessionLocal
from backend.db.models import ExpectedDamage


def yield_expected_damages(expected_fname):
    """Generate ExpectedDamage database model objects, reading through
    potentially very large parquet files in row batches.
    """
    df = pandas.read_parquet(expected_fname)

    batch_df = parse_exp_df(df)
    for row in batch_df.itertuples():
        if row.epoch == 1980:
            epoch = 2010
        else:
            epoch = row.epoch
        yield ExpectedDamage(
            feature_id=row.uid,
            hazard=row.hazard,
            rcp=row.rcp,
            epoch=epoch,
            protection_standard=row.protection_standard,
            ead_amin=row.ead_amin,
            ead_mean=row.ead_mean,
            ead_amax=row.ead_amax,
            eael_amin=row.eael_amin,
            eael_mean=row.eael_mean,
            eael_amax=row.eael_amax,
        )


def parse_exp_df(data):
    """Parse a fat df to long format"""
    # parquet files come in "fat" format, with string column names representing
    # damage type, defended status, and min/mean/max suffix
    data_cols = [c for c in data.columns if "EA" in c]

    # no protection_standard in list by default, append if in scheme
    id_vars = ["uid", "hazard", "rcp", "epoch"]
    if "protection_standard" in data_cols:
        id_vars.append("protection_standard")

    # melt to long format
    data = (
        data.melt(id_vars=id_vars, value_vars=data_cols)
        .query("value > 0")
        .reset_index(drop=True)
    )

    data.epoch = data.epoch.str.replace('hist', '1980').astype(int)

    # parse string key column for metadata
    meta = data.variable.str.extract(r"^([^_]+)_(\w+)_([^_]+)$")
    meta.columns = ["damage", "defended", "var"]

    # join metadata columns
    data = data.join(meta)
    if "protection_standard" not in data_cols:
        data["protection_standard"] = 0
    else:
        data.loc[data.defended == "undefended", "protection_standard"] = 0

    # pivot back up so we end with a row per uid, hazard etc. (see index columns below)
    # and columns for each damage type, each with min/mean/max
    data = (
        data.drop(columns=["variable", "var"])
        .pivot_table(
            index=["uid", "hazard", "rcp", "epoch", "protection_standard"],
            columns=["damage"],
            values="value",
        )
        .reset_index()
        .groupby(["uid", "hazard", "rcp", "epoch", "protection_standard"])
        .agg(['min', 'mean', 'max'])
        .fillna(0)
    )

    data.columns = [f"{var.lower()}_{stat}" for var, stat in data.columns]
    data = data.rename(
        columns={
            "ead_min": "ead_amin",
            "ead_mean": "ead_mean",
            "ead_max": "ead_amax",
            "eael_min": "eael_amin",
            "eael_mean": "eael_mean",
            "eael_max": "eael_amax",
        }
    )

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
        expected = snakemake.input.expected
        layer = snakemake.wildcards.layer
        output = snakemake.output
    except NameError:
        print("Expected to run from snakemake")
        exit()

    damage_exp = yield_expected_damages(expected)

    db: Session
    with SessionLocal() as db:
        for i, damage_exp in enumerate(tqdm(damage_exp, desc=f"{layer}_exp")):
            db.add(damage_exp)
            if i % 1000 == 0:
                db.commit()
        db.commit()

    with open(str(output), "w") as fh:
        fh.write(f"Loaded to database.\n\n")
        fh.write(f"From:\n{input}\n\n")
