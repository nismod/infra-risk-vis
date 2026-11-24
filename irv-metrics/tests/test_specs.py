"""Round-trip tests for irv_metrics.specs: upstream query/JSON path <-> canonical cube keys."""
import pytest

from irv_metrics.specs import DATASETS, MappingError, leaves


def load(dataset_id, kind, query, body):
    """Upstream response -> {(metric, dims): value}, as the loader flattens it."""
    dataset = DATASETS[dataset_id]
    ctx, assignment = dataset.apply(query)
    templates = dataset.profile if kind == "profile" else dataset.ladder
    out = {}
    for path, value, _ in leaves(body):
        for template in templates:
            found = template.match(list(path), ctx)
            if found is not None:
                metric, dims = dataset.canonical({**assignment, **found})
                out[(metric, tuple(sorted(dims.items())))] = value
                break
        else:
            raise AssertionError(f"unmapped {path}")
    return out


def test_risk_metric_splits_into_metric_and_return_period():
    ctx, assignment = DATASETS["risk/capital-stock"].apply({"risk_metric": "100", "sector": "total"})
    assert DATASETS["risk/capital-stock"].canonical(assignment) == ("loss", {"sector": "total", "return_period_years": "100"})
    ctx, assignment = DATASETS["risk/capital-stock"].apply({"risk_metric": "average_annual_loss_protected"})
    assert assignment == {"metric": "average_annual_loss"}


def test_population_group_dimension_depends_on_layer():
    dataset = DATASETS["risk/population"]
    facts = load("risk/population", "profile", {"population_layer": "wealth", "risk_metric": "50"},
                 {"adm_id": "KEN", "name": "Kenya", "groups": {"q1": {"value": 5, "unit": "people"}}})
    assert facts == {("exposure", (("return_period_years", "50"), ("wealth_group", "q1"))): 5}
    ctx, _ = dataset.apply({"population_layer": "demographics"})
    assert dataset.profile[0].render({"metric": "exposure", "demographic_group": "female"}, ctx) == ["groups", "female"]


def test_profile_leaves_with_embedded_members():
    facts = load("vulnerability/accessibility", "profile",
                 {"facility": "schools", "mode": "walking", "demographic_group": "total"},
                 {"adm_id": "X", "name": "X", "travel_time_avg": {"value": 12}, "pop_within_30min": {"value": 7}})
    base = (("demographic_group", "total"), ("facility_type", "schools"), ("travel_mode", "walking"))
    assert facts[("travel_time_avg", base)] == 12
    assert facts[("pop_within_threshold", tuple(sorted(base + (("travel_time_threshold_min", "30"),))))] == 7


def test_mapped_keys_round_trip():
    dataset = DATASETS["risk/social-facilities"]
    ctx, _ = dataset.apply({"infrastructure_type": "hospitals"})
    template = dataset.ladder[0]
    found = template.match(["thresholds", "moderate_plus", "100"], ctx)
    assert found == {"damage_threshold": "moderate", "return_period_years": "100", "metric": "count_rp_threshold"}
    assert template.render(found, ctx) == ["thresholds", "moderate_plus", "100"]


def test_metric_suffix_from_display_mode():
    dataset = DATASETS["hazard/tropical-cyclone-wind"]
    ctx, _ = dataset.apply({"display_mode": "percentage", "return_period_years": "10"})
    assert dataset.profile[1].match(["bands", "cat3"], ctx)["metric"] == "area_in_band_share"


def test_unknown_parameters_are_rejected():
    with pytest.raises(MappingError):
        DATASETS["exposure/roads"].apply({"colour": "red"})
    with pytest.raises(MappingError):
        DATASETS["exposure/population"].apply({"display_mode": "sideways"})


def test_specs_follow_the_config():
    """specs.py (the upstream mapping) uses the config's dataset ids and writes exactly its stored metrics,
    keyed by the same dimensions, so the scrape and CSV loaders give the same series."""
    from irv_metrics.config import load_catalogue
    from irv_metrics.specs import CURVES

    catalogue = load_catalogue()
    specs = {d.id: d for d in [*DATASETS.values(), CURVES]}
    assert specs.keys() == catalogue.datasets.keys()
    for dataset_id, spec in specs.items():
        stored = catalogue.datasets[dataset_id].stored_metrics
        assert {m: sorted(d) for m, d in spec.metrics.items()} == {m: sorted(s.dimensions) for m, s in stored.items()}, \
            dataset_id


def test_upstream_paths():
    """Only the curves are spelled differently upstream; every upstream path is distinct."""
    from irv_metrics.specs import CURVES

    assert CURVES.id == "adaptation-outcomes/concentration-curves" and CURVES.path == "concentration-curves"
    assert all(d.path == d.id for d in DATASETS.values())
    ctx, assignment = CURVES.apply({"hazard": "flood", "model": "jrc", "scenario": "baseline_protected"})
    assert CURVES.canonical(assignment) == ("cumulative_exposure_share", {"scenario": "baseline_protected"})
    with pytest.raises(MappingError):
        CURVES.apply({"hazard": "earthquake"})
