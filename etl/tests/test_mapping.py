"""layers.csv row -> catalogue item mapping."""

import pandas as pd
import pytest

import datapkg_mapping as dm
from conftest import AQUEDUCT_ITEMS, ETL


def test_parse_and_format_filter():
    f = dm.parse_filter("hazard=river; return_period=10|100;gcm=null")
    assert f == {"hazard": "river", "return_period": ["10", "100"], "gcm": None}
    assert dm.parse_filter(dm.format_filter(f)) == f
    assert dm.parse_filter("") == {}
    with pytest.raises(dm.MappingError):
        dm.parse_filter("hazard")


def test_matches_normalises_numbers_and_nulls():
    props = {"return_period": 100, "epoch": 2030, "gcm": None, "rcp": "4.5"}
    assert dm.matches(props, {"return_period": "100", "epoch": 2030.0, "rcp": 4.5})
    assert dm.matches(props, {"gcm": None, "missing": None})
    assert not dm.matches(props, {"gcm": "None"})
    assert dm.matches(props, {"return_period": ["10", "100"]})
    assert not dm.matches(props, {"return_period": 10})


def test_real_aqueduct_mapping_builds_filters():
    mapping = dm.load_mapping("aqueduct", ETL)
    layers = dm.read_layers("aqueduct", ETL)
    coastal = dm.layer_row(layers, "inuncoast_historical_wtsub_hist_rp0001_5")
    assert dm.row_filters(coastal, mapping) == {
        "hazard": "coastal", "rcp": "baseline", "epoch": "baseline", "gcm": None, "return_period": 1,
    }
    river = dm.layer_row(layers, "inunriver_rcp4p5_00000NorESM1-M_2030_rp00100")
    assert dm.row_filters(river, mapping)["hazard"] == "river"


def test_every_wired_pipeline_maps_every_row():
    for pipeline in ("aqueduct", "storm", "ghsl_pop", "jrc_flood", "landslide_arup"):
        mapping = dm.load_mapping(pipeline, ETL)
        layers = dm.read_layers(pipeline, ETL)
        specs = [dm.item_spec(r, mapping) for r in layers.to_dict(orient="records")]
        assert specs and all(s["dataset"] == mapping["dataset"] for s in specs)
        # mapped columns identify rows uniquely, so filters are distinct
        assert len({dm.format_filter(s["filters"]) for s in specs}) == len(specs), pipeline


def test_item_column_and_row_filter_override():
    mapping = {"dataset": "d", "properties": {"rp": "return_period"}, "values": {},
               "filters": {"hazard": "river"}}
    assert dm.item_spec({"rp": 10, "datapkg_item": "x"}, mapping) == {"dataset": "d", "item": "x"}
    assert dm.item_spec({"rp": 10, "datapkg_dataset": "e"}, mapping)["dataset"] == "e"
    spec = dm.item_spec({"rp": 10, "datapkg_filter": "hazard=coastal;slr=95"}, mapping)
    assert spec["filters"] == {"hazard": "coastal", "return_period": 10, "slr": "95"}


def test_resolve_against_catalogue(catalogue):
    cat = dm.open_catalog(str(catalogue), release="2026.1")
    mapping = dm.load_mapping("aqueduct", ETL)
    layers = dm.read_layers("aqueduct", ETL)
    for item_id, _, filename in AQUEDUCT_ITEMS:
        if filename:
            row = dm.layer_row(layers, filename[:-4])
            assert dm.resolve(cat, dm.item_spec(row, mapping)).id == item_id
    assert dm.resolve(cat, {"dataset": "aqueduct_flood", "item": item_id}).id == item_id
    with pytest.raises(dm.MappingError, match="matched 3 items"):
        dm.resolve(cat, {"dataset": "aqueduct_flood", "filters": {"hazard": "river"}})
    with pytest.raises(dm.MappingError, match="matched 0 items"):
        dm.resolve(cat, {"dataset": "aqueduct_flood", "filters": {"hazard": "pluvial"}})


def test_metadata_from_collection_keeps_local_fields():
    metadata = {"domain": "aqueduct", "name": "Local name", "group": "Hazard",
                "description": "old", "license": "old", "keys": ["hazard", "rp"]}
    collection = {"license": "CC-BY-4.0", "description": "Synthetic",
                  "sci:citation": "Ward et al. (2020)", "sci:doi": "10.1/x", "title": "Other"}
    out = dm.metadata_from_collection(metadata, collection)
    assert out == {**metadata, "license": "CC-BY-4.0",
                   "description": "Ward et al. (2020) https://doi.org/10.1/x"}
    assert dm.metadata_from_collection(metadata, collection, "description")["description"] == "Synthetic"
    # metadata.json is inserted as a raster_tile_sources row: no new keys
    assert set(out) == set(metadata)


def test_client_filters_accept_number_forms_and_absent():
    f = dm.client_filters({"return_period": 100.0, "gcm": None, "rcp": ["4.5", None]})
    assert f["return_period"] == ["100", "100.0"]
    assert f["gcm"] == [None] and f["rcp"] == [None, "4.5"]
