"""Mapping between the upstream GRI API and the canonical data cube.

The cube stores one value per (dataset, metric, {dimension: member}, adm_id).
The upstream API exposes the same values through three route kinds whose query
parameters and JSON layouts differ per dataset:

    series   /metrics/{dataset}/{iso3}/{level}?...       -> {unit, regions: [...]}
    profile  /metrics/{dataset}/{iso3}/profile/{adm_id}?  -> nested {key: {value, unit}}
    ladder   /metrics/{dataset}/{iso3}/ladder/{adm_id}?   -> {unit, return_periods: {...}}

Dataset ids are the cube's (contract C2) ids. Where an upstream route spells a
dataset differently, `Dataset.upstream` holds the upstream path (only the
concentration curves: `concentration-curves`); the scrape loader and the legacy
routes use `Dataset.path` at the upstream edge and the id everywhere else.

Each dataset is described once, declaratively, and the description is used in
both directions:

  * load (upstream -> canonical): query parameters give a partial assignment,
    each leaf's JSON path is matched against a template that supplies the rest;
  * serve (canonical -> upstream): the same parameters select observations and
    the same templates render each observation back to its JSON path.

Parameter specs:
    "dim"                      the parameter is the canonical dimension "dim"
    {"value": {...}, ...}      value -> assignment
    callable(value, ctx)       returns an assignment
Assignments are {dimension: member} plus optional "metric". Keys starting with
"$" are context variables (they steer templates; they are never stored).

Templates are dot-separated JSON paths. A segment may contain {var} or
{var:a|b} (restricted values); var is "metric", a dimension, or "$name" (a
dimension whose name comes from the context). Fixed assignments may use
"{$name}" placeholders.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import re
from typing import Callable

RP = "return_period_years"
FLOOD = "flood"

URBANISATION_CODES = {
    "urban_centre": "30", "dense_urban": "23", "semi_dense_urban_area": "22", "suburban_area": "21",
    "rural_settlement": "13", "low_density_rural": "12", "remote_area": "11",
}


class MappingError(ValueError):
    """A request or response that the mapping cannot interpret."""


class Template:
    _var = re.compile(r"\{(\$?\w+)(?::([\w|]+))?\}")

    def __init__(self, path: str, fixed: dict | None = None, maps: dict | None = None):
        self.path = path
        self.fixed = fixed or {}
        self.maps = maps or {}                       # var -> {api key: canonical member}
        self.unmaps = {v: {c: a for a, c in m.items()} for v, m in self.maps.items()}
        self.segments = []                           # (regex, [(var, domain)], format string)
        for segment in path.split("."):
            variables, regex, fmt, pos = [], "", "", 0
            for match in self._var.finditer(segment):
                literal = segment[pos:match.start()]
                regex += re.escape(literal)
                fmt += literal.replace("{", "{{").replace("}", "}}")
                name, domain = match.group(1), match.group(2)
                regex += f"({domain})" if domain else "(.+)"
                fmt += "{%d}" % len(variables)
                variables.append((name, set(domain.split("|")) if domain else None))
                pos = match.end()
            regex += re.escape(segment[pos:])
            fmt += segment[pos:].replace("{", "{{").replace("}", "}}")
            self.segments.append((re.compile(regex), variables, fmt))

    @staticmethod
    def resolve(name: str, ctx: dict) -> str:
        return ctx[name] if name.startswith("$") else name

    def free(self, ctx: dict) -> set[str]:
        return {self.resolve(name, ctx) for _, variables, _ in self.segments for name, _ in variables}

    def fixed_for(self, ctx: dict) -> dict:
        return {k: v.format(**{n[1:]: x for n, x in ctx.items()}) if "{" in v else v for k, v in self.fixed.items()}

    def match(self, keys: list[str], ctx: dict) -> dict | None:
        if len(keys) != len(self.segments):
            return None
        out = {}
        for key, (regex, variables, _) in zip(keys, self.segments):
            found = regex.fullmatch(str(key))
            if not found:
                return None
            for (name, _), raw in zip(variables, found.groups()):
                if name in self.maps:
                    if raw not in self.maps[name]:
                        return None
                    raw = self.maps[name][raw]
                out[self.resolve(name, ctx)] = raw
        return {**out, **self.fixed_for(ctx)}

    def render(self, fact: dict, ctx: dict) -> list[str] | None:
        keys = []
        for _, variables, fmt in self.segments:
            values = []
            for name, domain in variables:
                value = fact.get(self.resolve(name, ctx))
                if value is None:
                    return None
                value = self.unmaps.get(name, {}).get(value, value)
                if domain and value not in domain:
                    return None
                values.append(value)
            keys.append(fmt.format(*values))
        return keys


T = Template


@dataclass
class Dataset:
    id: str                                         # the cube's (C2) dataset id
    metrics: dict[str, tuple[str, ...]]             # metric -> dimensions it is keyed by
    params: dict[str, object] = field(default_factory=dict)
    series: dict = field(default_factory=dict)       # fixed assignment for series routes
    profile: list[Template] = field(default_factory=list)
    ladder: list[Template] = field(default_factory=list)
    echo: dict[str, tuple[str, ...]] = field(default_factory=dict)   # route kind -> echoed params
    nulls: tuple[str, ...] = ()                      # profile keys published as null
    ctx: dict = field(default_factory=dict)          # default context variables
    # display_mode=percentage: value as % of the sum of these assignments (derived, not stored).
    percent_of: tuple[dict, ...] = ()
    upstream: str | None = None                      # the upstream route's dataset path, where it differs from id

    @property
    def path(self) -> str:
        """The dataset's segment(s) in upstream /metrics/... routes."""
        return self.upstream or self.id

    def apply(self, query: dict[str, str]) -> tuple[dict, dict]:
        """Query parameters -> (context, partial canonical assignment)."""
        ctx, assignment, later = dict(self.ctx), {}, []
        for key, value in query.items():
            spec = self.params.get(key)
            if spec is None:
                raise MappingError(f"{self.id}: unknown parameter {key}")
            if callable(spec):
                later.append((spec, value))
                continue
            if isinstance(spec, str):
                part = {spec: value}
            elif value in spec:
                part = spec[value]
            else:
                raise MappingError(f"{self.id}: {key}={value} is not a known value")
            for k, v in part.items():
                (ctx if k.startswith("$") else assignment)[k] = v
        for spec, value in later:
            for k, v in spec(value, ctx).items():
                (ctx if k.startswith("$") else assignment)[k] = v
        return ctx, assignment

    def canonical(self, assignment: dict) -> tuple[str, dict] | None:
        """Assignment -> (metric, dims), dropping dimensions the metric is not keyed by."""
        metric = assignment.get("metric")
        if metric is None and len(self.metrics) == 1:
            metric = next(iter(self.metrics))
        if metric not in self.metrics:
            return None
        return metric, {d: str(assignment[d]) for d in self.metrics[metric] if d in assignment}

    def matches(self, metric: str, dims: dict, target: dict) -> bool:
        """Does an observation satisfy a target assignment (ignoring dims the metric lacks)?"""
        if "metric" in target and target["metric"] != metric:
            return False
        keyed = self.metrics[metric]
        return all(dims.get(k) == str(v) for k, v in target.items() if k != "metric" and k in keyed)


def risk_metric(rp_metric: str, annual_metric: str, annual_key: str):
    def apply(value, ctx):
        if value == annual_key:
            return {"metric": annual_metric}
        return {"metric": rp_metric, RP: value}
    return apply


def population_group(value, ctx):
    return {ctx["$group"]: value}


def travel_metric(value, ctx):
    found = re.fullmatch(r"pop_within_(\d+)min", value)
    if found:
        return {"metric": "pop_within_threshold", "travel_time_threshold_min": found.group(1)}
    return {"metric": value}


def statistic(value, ctx):
    if value == "avoided":
        return {"adaptation_state": "avoided"}
    return {"cost_bound": value}


PERCENT = {"absolute": {}, "percentage": {"$percent": "1"}}
NBS = {
    "opportunity_area": (), "biodiversity_benefit": (), "carbon_benefit": (),
    "implementation_cost": ("implementation_approach",),
}
NBS_PROFILE = [
    T("scoreboard.{metric:opportunity_area|biodiversity_benefit|carbon_benefit}"),
    T("scoreboard.{implementation_approach}", {"metric": "implementation_cost"}),
]
NBS_PARAMS = {"metric": "metric", "implementation_approach": "implementation_approach"}


def outcome(id: str, cost_metric: str, scenario_dims: tuple[str, ...], cost_bounds: bool) -> Dataset:
    s = scenario_dims
    metrics = {
        "average_annual_loss": s + ("sector", "adaptation_state"),
        "average_annual_flood_exposure": s + ("wealth_group", "adaptation_state"),
        "concentration_index": s + ("adaptation_state",),
        cost_metric: s + (("cost_bound",) if cost_bounds else ()),
    }
    profile = [
        T("economic.{sector}.{adaptation_state}", {"metric": "average_annual_loss"}),
        T("social.{wealth_group}.{adaptation_state}", {"metric": "average_annual_flood_exposure"}),
        T("concentration_index.{adaptation_state}", {"metric": "concentration_index"}),
    ]
    if cost_bounds:
        profile += [T(f"cost_{key}", {"metric": cost_metric, "cost_bound": bound})
                    for key, bound in (("quantity", "estimate"), ("lower", "lower"), ("upper", "upper"))]
    else:
        profile.append(T("cost_quantity", {"metric": cost_metric}))
    return Dataset(
        id, metrics,
        params={"metric": "metric", "statistic": statistic, "sector": "sector", "wealth_group": "wealth_group",
                **{d: d for d in s}},
        profile=profile,
        nulls=() if cost_bounds else ("cost_lower", "cost_upper"),
    )


DATASETS: dict[str, Dataset] = {d.id: d for d in [
    Dataset(
        "hazard/river-flooding",
        {m: (RP,) for m in ("flooded_area", "area_flooded_percentage", "mean_flood_depth", "p90_flood_depth")},
        params={"metric": "metric", RP: RP},
        profile=[T("metrics.{metric}")],
        ladder=[T("return_periods.{%s}" % RP)],
        echo={"profile": (RP,)},
    ),
    Dataset(
        "hazard/tropical-cyclone-wind",
        {"area_exceeding": (RP, "wind_threshold"), "area_in_band": (RP, "wind_band"),
         "area_exceeding_share": (RP, "wind_threshold"), "area_in_band_share": (RP, "wind_band")},
        params={RP: RP, "wind_threshold": "wind_threshold",
                "display_mode": {"area": {"$share": ""}, "percentage": {"$share": "_share"}}},
        ctx={"$share": ""},
        series={"metric": "area_exceeding{share}"},
        profile=[T("thresholds.{wind_threshold}", {"metric": "area_exceeding{share}"}),
                 T("bands.{wind_band}", {"metric": "area_in_band{share}"})],
        ladder=[T("return_periods.{%s}" % RP, {"metric": "area_exceeding{share}"})],
        echo={"profile": (RP, "display_mode")},
    ),
    Dataset(
        "exposure/population", {"population": ("demographic_group",)},
        params={"demographic_group": "demographic_group", "display_mode": PERCENT},
        profile=[T("groups.{demographic_group}")],
        percent_of=({"demographic_group": "total"},),
    ),
    Dataset("exposure/capital-stock", {"capital_stock": ("sector",)},
            params={"sector": "sector"}, profile=[T("sectors.{sector}")]),
    Dataset("exposure/roads", {"length": ("road_class",)},
            params={"road_class": "road_class"}, profile=[T("classes.{road_class}")]),
    Dataset("exposure/rail", {"length": ()}, profile=[T("length")]),
    Dataset("exposure/power", {"length": ()}, profile=[T("length")]),
    *[Dataset(
        f"exposure/{layer}", {"count": (), "per_100k": ("demographic_group",)},
        params={"metric": "metric", "demographic_group": "demographic_group"},
        profile=[T("{metric:count|per_100k}")],
    ) for layer in ("healthcare-facilities", "educational-facilities")],
    Dataset(
        "vulnerability/relative-wealth-index", {"rwi": ("rwi_measure",)},
        params={"rwi_measure": "rwi_measure"},
        profile=[T("{rwi_measure:average|population_weighted}")],
        echo={"series": ("rwi_measure",)},
    ),
    Dataset(
        "vulnerability/wealth-distribution", {"population": ("wealth_group",)},
        params={"wealth_group": "wealth_group", "display_mode": PERCENT},
        profile=[T("groups.{wealth_group}")],
        percent_of=tuple({"wealth_group": q} for q in ("q1", "q2", "q3", "q4", "q5")),
    ),
    Dataset(
        "vulnerability/accessibility",
        {m: ("facility_type", "travel_mode", "demographic_group") + extra
         for m, extra in (("travel_time_avg", ()), ("pop_within_threshold", ("travel_time_threshold_min",)))},
        params={"facility": "facility_type", "mode": "travel_mode", "demographic_group": "demographic_group",
                "metric": travel_metric},
        profile=[T("travel_time_avg", {"metric": "travel_time_avg"}),
                 T("pop_within_{travel_time_threshold_min}min", {"metric": "pop_within_threshold"})],
    ),
    Dataset(
        "risk/population",
        {"average_annual_exposure": ("demographic_group", "wealth_group"),
         "exposure": ("demographic_group", "wealth_group", RP)},
        params={"population_layer": {"demographics": {"$group": "demographic_group"}, "wealth": {"$group": "wealth_group"}},
                "population_group": population_group,
                "risk_metric": risk_metric("exposure", "average_annual_exposure", "average_annual_exposure_protected")},
        ctx={"$group": "demographic_group"},
        profile=[T("groups.{$group}")],
        ladder=[T("return_periods.{%s}" % RP, {"metric": "exposure"})],
    ),
    Dataset(
        "risk/capital-stock", {"average_annual_loss": ("sector",), "loss": ("sector", RP)},
        params={"sector": "sector",
                "risk_metric": risk_metric("loss", "average_annual_loss", "average_annual_loss_protected")},
        profile=[T("sectors.{sector}")],
        ladder=[T("return_periods.{%s}" % RP, {"metric": "loss"})],
    ),
    Dataset(
        "risk/direct-damage", {"expected_annual_damage": ("hazard", "infrastructure_type", "road_class")},
        params={"hazard": "hazard", "infrastructure_type": "infrastructure_type", "asset_class": "road_class"},
        profile=[T("classes.{road_class}")],
    ),
    Dataset(
        "risk/social-facilities",
        {"ea_damages": ("facility_type",), "count_rp_threshold": ("facility_type", RP, "damage_threshold")},
        params={"infrastructure_type": "facility_type", "metric": "metric", RP: RP, "damage_threshold": "damage_threshold"},
        profile=[T("ea_damages", {"metric": "ea_damages"}),
                 T("thresholds.{damage_threshold}", {"metric": "count_rp_threshold"},
                   maps={"damage_threshold": {"moderate_plus": "moderate", "extensive_plus": "extensive", "complete": "complete"}})],
        ladder=[T("thresholds.{damage_threshold}.{%s}" % RP, {"metric": "count_rp_threshold"},
                  maps={"damage_threshold": {"moderate_plus": "moderate", "extensive_plus": "extensive", "complete": "complete"}})],
        echo={"profile": (RP,)},
    ),
    Dataset(
        "risk/social-accessibility",
        {m: ("facility_type", "travel_mode", "demographic_group", "disruption_type") for m in ("ea_people", "ea_person_min")},
        params={"infrastructure_type": "facility_type", "mode": "travel_mode", "demographic_group": "demographic_group",
                "disruption_type": "disruption_type",
                "measure": {"people": {"metric": "ea_people"}, "travel_time_increase": {"metric": "ea_person_min"}}},
        profile=[T("ea_people_{disruption_type}", {"metric": "ea_people"}),
                 T("ea_person_min_{disruption_type}", {"metric": "ea_person_min"})],
    ),
    Dataset(
        "adaptation-potential/slope-vegetation", {m: d + ("land_use",) for m, d in NBS.items()},
        params={**NBS_PARAMS, "land_use": "land_use"},
        profile=NBS_PROFILE + [T("classes.{land_use}")],
        echo={"profile": ("land_use",)},
    ),
    Dataset(
        "adaptation-potential/mangroves", {m: d + ("shoreline_condition",) for m, d in NBS.items()},
        params={**NBS_PARAMS, "shoreline_condition": "shoreline_condition"},
        profile=NBS_PROFILE + [T("classes.{shoreline_condition}")],
        echo={"profile": ("shoreline_condition",)},
    ),
    Dataset("adaptation-potential/river-catchment-restoration", dict(NBS), params=NBS_PARAMS, profile=NBS_PROFILE),
    Dataset("adaptation-potential/existing-flood-protection", {"protection_standard": ()},
            profile=[T("protection_standard_mode")]),
    Dataset("adaptation-potential/river-network-context", {"river_length": ("urbanisation_class",)},
            params={"urbanisation_class": "urbanisation_class"},
            profile=[T("urbanisation_class.{urbanisation_class}")]),
    outcome("adaptation-outcomes/dry-proofing", "area_dry_proofed", (), cost_bounds=False),
    outcome("adaptation-outcomes/relocation", "capstock_relocated", ("urbanisation_threshold",), cost_bounds=False),
    outcome("adaptation-outcomes/flood-protection", "adaptation_cost",
            ("urbanisation_threshold", "design_return_period_years"), cost_bounds=True),
]}

# Concentration curves: national Lorenz-type curves of flood exposure by wealth. Upstream serves them at
# /metrics/concentration-curves/{iso3}; the cube (C2) id is adaptation-outcomes/concentration-curves. The
# upstream `hazard=flood` and `model=jrc` are the same for every curve, so they are not dimensions.
CURVES = Dataset(
    "adaptation-outcomes/concentration-curves",
    {"cumulative_exposure_share": ("scenario", "urbanisation_threshold", "design_return_period_years")},
    upstream="concentration-curves",
    params={"hazard": {"flood": {}}, "model": {"jrc": {}}, "scenario": "scenario",
            "urbanisation_threshold_code": lambda code, ctx: {
                "urbanisation_threshold": {c: u for u, c in URBANISATION_CODES.items()}[code]},
            "design_return_period_years": "design_return_period_years"},
)
CURVE_ECHO = ("hazard", "model", "scenario")

# Keys in profile/ladder bodies that identify the place or echo the request.
STRUCTURAL = {"adm_id", "name", "unit"}


def leaves(body, path=()):
    """Yield (path, value, unit) for every measured value in a profile/ladder body."""
    for key, item in body.items():
        if not path and key in STRUCTURAL:
            continue
        here = path + (str(key),)
        if isinstance(item, dict) and "value" in item:
            yield here, item["value"], item.get("unit")
        elif isinstance(item, dict):
            yield from leaves(item, here)
        elif isinstance(item, (int, float)) and not isinstance(item, bool):
            yield here, item, None
        elif item is None:
            yield here, None, None


def fill(assignment: dict, ctx: dict) -> dict:
    """Substitute {name} context placeholders ($name in ctx) in a fixed assignment."""
    names = {k[1:]: v for k, v in ctx.items()}
    return {k: v.format(**names) if isinstance(v, str) and "{" in v else v for k, v in assignment.items()}
