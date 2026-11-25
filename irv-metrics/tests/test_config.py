"""The dataset config under config/gri loads, validates against contract C2 and cross-checks."""
from __future__ import annotations

import shutil

import pytest
import yaml

from irv_metrics import config
from irv_metrics.config import ConfigError, load_catalogue

EXAMPLES = config.PACKAGE.parent / "datasets" / "examples"


@pytest.fixture(scope="module")
def catalogue():
    return load_catalogue()


def test_gri_config_loads(catalogue):
    assert list(catalogue.sections) == ["hazard", "exposure", "vulnerability", "risk",
                                        "adaptation-potential", "adaptation-outcomes"]
    assert len(catalogue.datasets) == 26
    assert catalogue.boundary_sets["cgaz-6.0.0"].kind == "admin"
    assert catalogue.boundary_sets["cgaz-6.0.0"].level_numbers == [0, 1, 2]


def test_every_dataset_reads_published_csvs(catalogue):
    for dataset in catalogue.datasets.values():
        assert dataset.sources, dataset.id
        for source in dataset.sources:
            assert "boundary_set" not in source, f"{dataset.id}: GRI sources are on each country's admin set"


def test_datasets_in_display_order(catalogue):
    ids = list(catalogue.datasets)
    assert ids[:2] == ["hazard/river-flooding", "hazard/tropical-cyclone-wind"]
    assert ids[-1] == "adaptation-outcomes/concentration-curves"
    sections = [catalogue.datasets[i].section for i in ids]
    assert sections == sorted(sections, key=list(catalogue.sections).index)


def test_typed_access(catalogue):
    population = catalogue.datasets["exposure/population"]
    assert list(population.stored_metrics) == ["population"]
    assert population.metrics["population_share"].derived
    assert population.metric_dimensions("population_share") == ("demographic_group",)
    assert population.dimensions["demographic_group"].default == "total"
    assert catalogue.members("hazard/river-flooding", "return_period_years") == ("10", "20", "50", "75", "100", "200", "500")
    assert catalogue.dimensions["urbanisation_threshold"].by_code("30") == "urban_centre"
    assert catalogue.member_label("exposure/population", "demographic_group", "infant") == "Children Under 5"
    curves = catalogue.datasets["adaptation-outcomes/concentration-curves"]
    assert curves.kind == "curve" and curves.x["id"] == "cumulative_population_share"


@pytest.mark.parametrize("name", ["gri", "standalone-jam"])
def test_examples_load(name):
    """The loader reads the C2 worked examples too (overlay and group sets, a range dimension, wide sources)."""
    load_catalogue(EXAMPLES / name)


def test_schema_dir_argument_and_env(tmp_path, monkeypatch):
    copy = tmp_path / "schema"
    shutil.copytree(config.SCHEMA_DIR, copy)
    load_catalogue(schema_dir=copy)
    monkeypatch.setenv(config.SCHEMA_ENV, str(copy))
    assert config.schema_dir() == copy
    load_catalogue()
    with pytest.raises(FileNotFoundError, match="IRV_METRICS_SCHEMA_DIR"):
        load_catalogue(schema_dir=tmp_path / "missing")


def _broken(tmp_path, dataset: str, edit) -> ConfigError:
    root = tmp_path / "gri"
    shutil.copytree(config.CONFIG_DIR, root)
    path = root / f"{dataset}.yaml"
    document = yaml.safe_load(path.read_text())
    edit(document)
    path.write_text(yaml.safe_dump(document, sort_keys=False))
    with pytest.raises(ConfigError) as caught:
        load_catalogue(root)
    return caught.value


def test_schema_errors_are_reported(tmp_path):
    error = _broken(tmp_path, "hazard/river-flooding", lambda d: d["dimensions"].update(level={}))
    assert any("river-flooding.yaml: dimensions" in p for p in error.problems), error.problems


def test_cross_check_errors_are_reported(tmp_path):
    def edit(d):
        d["metrics"]["flooded_area"]["unit"] = "furlongs"
        d["dimensions"]["return_period_years"]["default"] = "1000"
        d["display"]["map"]["dims"]["return_period_years"] = "3"
        d["id"] = "hazard/other"
    error = _broken(tmp_path, "hazard/river-flooding", edit)
    text = "\n".join(error.problems)
    assert "unknown unit furlongs" in text
    assert "default 1000 is not one of the dataset's members" in text
    assert "display.map: return_period_years=3 is not a member" in text
    assert "does not match the file path" in text


def test_repeated_vocabulary_id_is_an_error(tmp_path):
    root = tmp_path / "gri"
    shutil.copytree(config.CONFIG_DIR, root)
    (root / "vocabulary-extra.yaml").write_text("units:\n  km: {label: km, quantity: length}\n")
    with pytest.raises(ConfigError, match="units.km is already defined"):
        load_catalogue(root)


class Recorder:
    """A cursor that records what is written: table -> rows."""

    def __init__(self):
        self.rows: dict[str, list[tuple]] = {}

    def execute(self, sql, params=None):
        if sql.startswith("INSERT"):
            self.rows.setdefault(sql.split()[2], []).append(params)

    def executemany(self, sql, rows):
        self.rows.setdefault(sql.split()[2], []).extend(rows)

    def fetchall(self):
        return []


def test_seed_writes_the_catalogue_under_its_ids(catalogue):
    """seed_catalogue writes every dataset of the config, under its C2 id, with sections from the config."""
    from irv_metrics.loaders import common

    cursor = Recorder()
    common.seed_catalogue(cursor, catalogue)
    rows = cursor.rows
    assert [r[0] for r in rows["section"]] == list(catalogue.sections)
    assert [r[0] for r in rows["dataset"]] == list(catalogue.datasets)            # display order, C2 ids
    assert "adaptation-outcomes/concentration-curves" in [r[0] for r in rows["dataset"]]
    assert {r[0] for r in rows["unit"]} == set(catalogue.units)
    members = {(r[0], r[1]): r for r in rows["dimension_member"]}
    assert members[("urbanisation_threshold", "urban_centre")][7] == "30"         # code
    assert members[("return_period_years", "1000")][6] == 1000                    # value
    assert members[("return_period_years", "10")][5] == 0                         # display order
    metrics = {(r[0], r[1]): r for r in rows["metric"]}
    assert metrics[("hazard/river-flooding", "flooded_area")][3] == "km2"
    assert ("exposure/population", "population_share") not in metrics            # share_of: computed by the API
    assert ("hazard/tropical-cyclone-wind", "area_in_band") in metrics            # stored, computed at load
    defaults = {(r[0], r[1]): r[2] for r in rows["dataset_dimension"]}
    assert defaults[("exposure/population", "demographic_group")] == "total"
    assert not any(r[1] == "hazard" for r in rows["dataset_dimension"] if r[0].endswith("concentration-curves"))
    assert [r[0] for r in rows["boundary_set"]] == ["cgaz-6.0.0"]
    assert [r[1:] for r in rows["boundary_set_level"]] == [(0, "Country"), (1, "Admin 1"), (2, "Admin 2")]


def test_seed_takes_overlay_sets_after_their_admin_set():
    """Boundary sets are registered admin first, whatever their order in the vocabulary (countries_from)."""
    from irv_metrics.loaders import common

    cursor = Recorder()
    common.seed_catalogue(cursor, load_catalogue(EXAMPLES / "gri"))
    sets = [(r[0], r[1]) for r in cursor.rows["boundary_set"]]
    kinds = [kind for _, kind in sets]
    assert kinds == sorted(kinds, key=lambda k: k != "admin") and "overlay" in kinds


def test_check_series(catalogue):
    from irv_metrics.loaders.common import check_series

    flood = "hazard/river-flooding"
    assert check_series(catalogue, [(flood, "flooded_area", (("return_period_years", "10"),))]) == []
    assert check_series(catalogue, [
        ("concentration-curves", "cumulative_exposure_share", ()),
        (flood, "flooded_area", (("return_period_years", "1000"),)),
        (flood, "flooded_area", (("hazard", "river_flood"),)),
        ("exposure/population", "population_share", ()),
    ]) == ["concentration-curves: not in the config",
           "exposure/population: population_share is not a stored metric",
           f"{flood}: flooded_area is not keyed by hazard (=river_flood)",
           f"{flood}: return_period_years=1000 is not one of the dataset's members"]
