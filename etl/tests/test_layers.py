"""layers_from_catalogue.py and metadata_from_catalogue.py."""

import json

import pandas as pd

import datapkg_mapping as dm
import layers_from_catalogue as lfc
import metadata_from_catalogue as mfc
from conftest import AQUEDUCT_ITEMS


def test_refresh_fills_item_ids_and_keeps_columns(catalogue, etl_root):
    before = dm.read_layers("aqueduct", etl_root)
    rc = lfc.main(["aqueduct", "--root", str(etl_root), "--base", str(catalogue),
                   "--release", "2026.1"])
    assert rc == 0
    after = dm.read_layers("aqueduct", etl_root)
    assert list(after.columns) == list(before.columns) + ["datapkg_item"]
    pd.testing.assert_frame_equal(after[before.columns], before)
    expected = {f: i for i, _, f in AQUEDUCT_ITEMS if f}
    assert dict(zip(after.filename, after.datapkg_item)) == expected
    # idempotent, and --check passes once refreshed
    assert lfc.main(["aqueduct", "--root", str(etl_root), "--base", str(catalogue),
                     "--release", "2026.1", "--check"]) == 0


def test_add_missing_appends_rows_from_item_properties(catalogue, etl_root):
    cat = dm.open_catalog(str(catalogue), release="2026.1")
    layers, report = lfc.generate("aqueduct", cat, etl_root, add_missing=True)
    assert report["new"] == ["inunriver_rcp8p5_noresm1-m_2080_rp01000.tif"]
    new = layers[layers.filename == report["new"][0]].iloc[0]
    assert new.hazard == "fluvial"  # reverse of the river value map
    assert (new.rp, new.rcp, new.epoch, new.gcm) == (1000, "8.5", 2080, "NorESM1-M")
    assert new.key == "inunriver_rcp8p5_noresm1-m_2080_rp01000"
    assert new["var"] == "depth_fluvial"
    assert new.datapkg_item == "inunriver_rcp8p5_noresm1-m_2080_rp01000"


def test_generate_from_scratch(catalogue, etl_root):
    (etl_root / "pipelines/storm/layers.csv").unlink()
    cat = dm.open_catalog(str(catalogue), release="2026.1")
    layers, report = lfc.generate("storm", cat, etl_root)
    assert sorted(layers.filename) == ["cmcc_2050_rp10.tif", "present_rp10.tif"]
    assert set(layers.columns) == {"filename", "rp", "gcm", "rcp", "epoch", "var", "datapkg_item"}
    assert set(layers["var"]) == {"ws"}


def test_check_reports_unmatched(catalogue, etl_root, capsys):
    path = etl_root / "pipelines/aqueduct/layers.csv"
    layers = pd.read_csv(path, keep_default_na=False)
    layers.loc[0, "rp"] = 7
    layers.to_csv(path, index=False)
    rc = lfc.main(["aqueduct", "--root", str(etl_root), "--base", str(catalogue), "--check"])
    assert rc == 1
    assert "unmatched row" in capsys.readouterr().err


def test_metadata_from_catalogue(catalogue, etl_root):
    path = etl_root / "pipelines/aqueduct/metadata.json"
    before = json.loads(path.read_text())
    assert mfc.main(["aqueduct", "--root", str(etl_root), "--base", str(catalogue),
                     "--check"]) == 1
    assert mfc.main(["aqueduct", "--root", str(etl_root), "--base", str(catalogue)]) == 0
    after = json.loads(path.read_text())
    for key in ("domain", "name", "group", "keys"):
        assert after[key] == before[key]
    assert after["license"] == "CC-BY-4.0"
    assert after["description"].startswith("Nobody (2026) Synthetic aqueduct.")
    assert mfc.main(["aqueduct", "--root", str(etl_root), "--base", str(catalogue),
                     "--check"]) == 0
