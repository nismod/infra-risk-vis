"""Load expected damages from a single file to database
"""
import pyarrow.parquet as pq

from sqlalchemy.orm import Session
from tqdm import tqdm

from backend.db.database import SessionLocal
from backend.db.models import ExpectedDamage


def yield_expected_damages(expected_fname):
    """Generate ExpectedDamage database model objects, reading through
    potentially very large parquet files in row batches.
    """
    pf = pq.ParquetFile(expected_fname)
    batch_size = 100
    batches = pf.iter_batches(batch_size)

    for batch in batches:
        batch_df = parse_exp_damage_batch(batch)
        for row in batch_df.itertuples():
            yield ExpectedDamage(
                feature_id=row.uid,
                hazard=row.hazard,
                rcp=row.rcp,
                epoch=row.epoch,
                protection_standard=0,
                ead_amin=row.ead_amin,
                ead_mean=row.ead_mean,
                ead_amax=row.ead_amax,
                eael_amin=0,
                eael_mean=0,
                eael_amax=0,
            )


def parse_exp_damage_batch(batch):
    """Parse a parquet (arrow) row batch to pandas"""
    data = batch.to_pandas()

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
