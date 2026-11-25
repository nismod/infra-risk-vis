"""Dataset config (contract C2): load, validate, cross-check, and read source files.

A config directory holds vocabulary files (`vocabulary*.yaml`, merged) and one
file per dataset at `{section}/{dataset}.yaml`:

    catalogue = load_catalogue(CONFIG_DIR)
    catalogue.datasets["hazard/river-flooding"].metrics["flooded_area"].unit   # 'km2'

Every file is validated against the C2 JSON Schemas, then cross-checked
(sections, units, dimensions, members, defaults, derived metrics, display
hints, boundary sets). All problems are collected and raised together as a
`ConfigError`.

`read_source` applies a dataset's `source` rules to the files under a data
directory and yields one result per value: `Placed`, `Skipped` or `Unplaced`.
It is the rule semantics of `datasets/examples/check.py`, made reusable for the
CSV loader.

The catalogue is the only definition of the cube's structure:
`loaders.common.seed_catalogue` writes it to the database under these ids.
"""
from __future__ import annotations

from collections.abc import Iterator
import csv
from dataclasses import dataclass, field
import json
import os
from pathlib import Path
import re

import jsonschema
import yaml

PACKAGE = Path(__file__).resolve().parent
# config/ and datasets/ sit next to the package in a source checkout (and an
# editable install); a built wheel bundles them as irv_metrics/_data/ (setup.py).
DATA_ROOT = PACKAGE / "_data" if (PACKAGE / "_data").is_dir() else PACKAGE.parent
CONFIG_DIR = DATA_ROOT / "config" / "gri"
# The schemas live with the frozen contract. A copy of the package elsewhere
# passes `schema_dir` or sets IRV_METRICS_SCHEMA_DIR.
SCHEMA_DIR = DATA_ROOT / "datasets" / "schema"
SCHEMA_ENV = "IRV_METRICS_SCHEMA_DIR"
DATASET_SCHEMA = "dataset.schema.json"
VOCABULARY_SCHEMA = "vocabulary.schema.json"

_TEMPLATE = re.compile(r"\{(\w+)(?:\|(\w+))?\}")


class ConfigError(ValueError):
    """The config does not validate or does not hang together. `problems` lists every one."""

    def __init__(self, problems: list[str]):
        self.problems = problems
        super().__init__(f"{len(problems)} config problem(s):\n" + "\n".join(f"  {p}" for p in problems))


# --- the catalogue -----------------------------------------------------------

@dataclass(frozen=True)
class Section:
    id: str
    title: str
    description: str | None = None


@dataclass(frozen=True)
class Unit:
    id: str
    label: str
    quantity: str


@dataclass(frozen=True)
class Member:
    id: str
    label: str
    short_label: str | None = None
    description: str | None = None
    value: float | None = None
    code: str | None = None


@dataclass(frozen=True)
class Dimension:
    id: str
    label: str
    members: dict[str, Member]                    # in display order
    description: str | None = None
    ordered: bool = False
    role: str | None = None

    def by_code(self, code: str) -> str:
        for member in self.members.values():
            if member.code == code:
                return member.id
        raise KeyError(f"{self.id}: no member with code {code!r}")


@dataclass(frozen=True)
class BoundarySet:
    id: str
    title: str
    kind: str
    levels: tuple[dict, ...]
    countries: dict | None = None
    attribution: str | None = None
    attribution_url: str | None = None

    @property
    def level_numbers(self) -> list[int]:
        return [entry.get("level", i) for i, entry in enumerate(self.levels)]


@dataclass(frozen=True)
class Metric:
    id: str
    label: str
    unit: str
    dimensions: tuple[str, ...] = ()
    aggregation: str = "none"
    weight: str | None = None
    description: str | None = None
    source_note: str | None = None
    derive: dict | None = None

    @property
    def derived(self) -> bool:
        return self.derive is not None


@dataclass(frozen=True)
class DatasetDimension:
    id: str
    members: tuple[str, ...]                      # the dataset's members, in display order
    default: str | None = None
    labels: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class Dataset:
    id: str
    title: str
    metrics: dict[str, Metric]                    # in display order
    dimensions: dict[str, DatasetDimension]
    sources: tuple[dict, ...] = ()                # C2 `source` entries, as written
    display: dict = field(default_factory=dict)
    kind: str = "values"
    x: dict | None = None
    description: str | None = None
    documentation: str | None = None
    source_note: str | None = None
    order: int | None = None
    path: Path | None = None

    @property
    def section(self) -> str:
        return self.id.split("/")[0]

    @property
    def stored_metrics(self) -> dict[str, Metric]:
        """Metrics whose values are loaded. `aggregate_of` metrics are stored too, but computed at load."""
        return {k: m for k, m in self.metrics.items() if not m.derived}

    def metric_dimensions(self, metric: str) -> tuple[str, ...]:
        """A metric's dimensions; a derived metric inherits them from the metric it derives from."""
        spec = self.metrics[metric]
        if spec.derive and "share_of" in spec.derive:
            return self.metrics[spec.derive["share_of"]].dimensions
        if spec.derive and "aggregate_of" in spec.derive:
            return self.metrics[spec.derive["aggregate_of"]].dimensions
        return spec.dimensions

    def boundary_sets(self, default_set: str | None = None) -> list[str]:
        """The dataset's sets in source order. A source without `boundary_set` is the country's admin set."""
        return [s.get("boundary_set", default_set) for s in self.sources]


@dataclass(frozen=True)
class Catalogue:
    sections: dict[str, Section]
    units: dict[str, Unit]
    dimensions: dict[str, Dimension]
    boundary_sets: dict[str, BoundarySet]
    countries: dict[str, dict]
    datasets: dict[str, Dataset]                  # by section order, then `order`, then id

    def members(self, dataset_id: str, dimension: str) -> tuple[str, ...]:
        return self.datasets[dataset_id].dimensions[dimension].members

    def member_label(self, dataset_id: str, dimension: str, member: str) -> str:
        override = self.datasets[dataset_id].dimensions[dimension].labels.get(member)
        return override or self.dimensions[dimension].members[member].label


# --- loading -----------------------------------------------------------------

def schema_dir(path: str | Path | None = None) -> Path:
    return Path(path or os.environ.get(SCHEMA_ENV) or SCHEMA_DIR)


def validators(path: str | Path | None = None) -> tuple[jsonschema.Validator, jsonschema.Validator]:
    """(dataset validator, vocabulary validator) from the C2 schema files."""
    root = schema_dir(path)
    out = []
    for name in (DATASET_SCHEMA, VOCABULARY_SCHEMA):
        file = root / name
        if not file.is_file():
            raise FileNotFoundError(f"C2 schema {file} not found; pass schema_dir or set {SCHEMA_ENV}")
        schema = json.loads(file.read_text())
        cls = jsonschema.validators.validator_for(schema)
        cls.check_schema(schema)
        out.append(cls(schema, format_checker=cls.FORMAT_CHECKER))
    return out[0], out[1]


def _read_yaml(path: Path, problems: list[str], name: str):
    try:
        return yaml.safe_load(path.read_text())
    except yaml.YAMLError as error:
        problems.append(f"{name}: {error}")
        return None


def _schema_problems(validator, document, name: str) -> list[str]:
    return [f"{name}: {'/'.join(map(str, e.absolute_path)) or '(root)'}: {e.message}"
            for e in sorted(validator.iter_errors(document), key=lambda e: list(map(str, e.absolute_path)))]


def load_catalogue(config_dir: str | Path = CONFIG_DIR, schema_dir: str | Path | None = None) -> Catalogue:
    """Load every file under `config_dir`, validate and cross-check it. Raises ConfigError listing all problems."""
    root = Path(config_dir)
    if not root.is_dir():
        raise ConfigError([f"{root}: not a directory"])
    dataset_validator, vocabulary_validator = validators(schema_dir)
    problems: list[str] = []
    files = sorted(root.rglob("*.yaml")) + sorted(root.rglob("*.yml"))
    vocab_files = [p for p in files if p.name.startswith("vocabulary")]
    dataset_files = [p for p in files if not p.name.startswith("vocabulary")]
    if not vocab_files:
        problems.append(f"{root}: no vocabulary*.yaml file")

    vocab: dict[str, dict] = {"sections": {}, "units": {}, "dimensions": {}, "boundary_sets": {}, "countries": {}}
    for path in vocab_files:
        name = str(path.relative_to(root))
        document = _read_yaml(path, problems, name)
        if document is None:
            continue
        found = _schema_problems(vocabulary_validator, document, name)
        problems += found
        if found:
            continue
        for part, entries in document.items():
            for key, value in entries.items():
                if key in vocab[part]:
                    problems.append(f"{name}: {part}.{key} is already defined in another vocabulary file")
                vocab[part][key] = value
    _check_vocabulary(vocab, problems)

    raw: dict[str, tuple[Path, dict]] = {}
    for path in dataset_files:
        name = str(path.relative_to(root))
        document = _read_yaml(path, problems, name)
        if document is None:
            continue
        found = _schema_problems(dataset_validator, document, name)
        problems += found
        if found:
            continue
        expected = path.relative_to(root).with_suffix("").as_posix()
        if document["id"] != expected:
            problems.append(f"{name}: id {document['id']} does not match the file path ({expected})")
        if document["id"] in raw:
            problems.append(f"{name}: dataset {document['id']} is defined twice")
        raw[document["id"]] = (path, document)
        _cross_check(name, document, vocab, problems)

    if problems:
        raise ConfigError(problems)
    return _build(vocab, raw)


def _members_of(dimension: dict) -> dict:
    if "range" in dimension:
        r = dimension["range"]
        return {str(i): {"label": str(i), "value": i} for i in range(r["from"], r["to"] + 1, r.get("step", 1))}
    return dimension["members"]


def _sources(document: dict) -> list[dict]:
    source = document.get("source")
    return [] if source is None else source if isinstance(source, list) else [source]


def _level_numbers(boundary_set: dict) -> list[int]:
    return [entry.get("level", i) for i, entry in enumerate(boundary_set["levels"])]


def _check_vocabulary(vocab: dict, problems: list[str]):
    sets = vocab["boundary_sets"]
    for set_id, spec in sets.items():
        numbers = _level_numbers(spec)
        if len(set(numbers)) != len(numbers) or numbers != sorted(numbers):
            problems.append(f"vocabulary: {set_id}: levels must be distinct and increasing")
        kind = spec.get("kind", "admin")
        if kind == "admin" and numbers[0] != 0:
            problems.append(f"vocabulary: {set_id}: an admin set starts at level 0, the country")
        if kind == "overlay" and sets.get(spec["countries"]["from"], {}).get("kind", "admin") != "admin":
            problems.append(f"vocabulary: {set_id}: countries.from must name an admin set")
        if kind == "overlay" and spec["countries"]["from"] not in sets:
            problems.append(f"vocabulary: {set_id}: unknown set {spec['countries']['from']} in countries.from")
    for iso3, spec in vocab["countries"].items():
        if "boundary_set" in spec and sets.get(spec["boundary_set"], {}).get("kind", "admin") != "admin":
            problems.append(f"vocabulary: countries.{iso3}: boundary_set must be an admin set")
        if "boundary_set" in spec and spec["boundary_set"] not in sets:
            problems.append(f"vocabulary: countries.{iso3}: unknown boundary set {spec['boundary_set']}")
    for dim_id, dim in vocab["dimensions"].items():
        codes = [m["code"] for m in _members_of(dim).values() if "code" in m]
        if len(codes) != len(set(codes)):
            problems.append(f"vocabulary: {dim_id}: member codes must be unique")


def _cross_check(name: str, d: dict, vocab: dict, problems: list[str]):
    """The references in one dataset file. The checks of examples/check.py, plus defaults and display hints."""
    def fail(message):
        problems.append(f"{name}: {message}")

    if d["id"].split("/")[0] not in vocab["sections"]:
        fail(f"unknown section in {d['id']}")
    used = d.get("dimensions", {})
    allowed = {}
    for dim, spec in used.items():
        if dim not in vocab["dimensions"]:
            fail(f"unknown dimension {dim}")
            continue
        known = _members_of(vocab["dimensions"][dim])
        for member in [*spec.get("members", []), *filter(None, [spec.get("default")]), *spec.get("labels", {})]:
            if member not in known:
                fail(f"{dim}: unknown member {member}")
        allowed[dim] = set(spec.get("members") or known)
        if "default" in spec and spec["default"] not in allowed[dim]:
            fail(f"{dim}: default {spec['default']} is not one of the dataset's members")
    metrics = d["metrics"]
    for metric, spec in metrics.items():
        if spec["unit"] not in vocab["units"]:
            fail(f"{metric}: unknown unit {spec['unit']}")
        for dim in spec.get("dimensions", []):
            if dim not in used:
                fail(f"{metric}: dimension {dim} is not listed under dimensions")
        if spec.get("aggregation") == "weighted_mean" and "weight" not in spec:
            fail(f"{metric}: weighted_mean needs a weight")
        if "weight" in spec and not re.fullmatch(r"[a-z][a-z0-9-]*/[a-z][a-z0-9-]*:[a-z][a-z0-9_]*", spec["weight"]):
            fail(f"{metric}: weight {spec['weight']} is not {{section}}/{{dataset}}:{{metric}}")
        derive = spec.get("derive", {})
        if "share_of" in derive:
            base = metrics.get(derive["share_of"])
            if base is None or "derive" in base:
                fail(f"{metric}: share_of must name a stored metric")
            else:
                for assignment in derive["total"]:
                    for dim, member in assignment.items():
                        if dim not in base.get("dimensions", []):
                            fail(f"{metric}: total names {dim}, which {derive['share_of']} is not keyed by")
                        elif dim in allowed and member not in allowed[dim]:
                            fail(f"{metric}: total names {dim}={member}, not a member")
        if "aggregate_of" in derive:
            base = metrics.get(derive["aggregate_of"])
            if base is None or "derive" in base:
                fail(f"{metric}: aggregate_of must name a stored metric")
            elif base.get("aggregation", "none") == "none":
                fail(f"{metric}: {derive['aggregate_of']} has aggregation none")
            elif base.get("aggregation") == "weighted_mean" and "weight" not in base:
                fail(f"{metric}: {derive['aggregate_of']} is a weighted mean without a weight")
            for target in derive["to"]:
                if vocab["boundary_sets"].get(target, {}).get("kind") != "group":
                    fail(f"{metric}: {target} is not a group set")
    if d.get("kind") == "curve" and d.get("x", {}).get("unit") and d["x"]["unit"] not in vocab["units"]:
        fail(f"x: unknown unit {d['x']['unit']}")

    display = d.get("display", {})
    for control in display.get("controls", []):
        if control != "metric" and control not in used:
            fail(f"display.controls: unknown dimension {control}")
    for chart in display.get("charts", []):
        for key in ("by", "split"):
            if chart.get(key, "area") != "area" and chart[key] not in used:
                fail(f"display.charts: unknown dimension {chart[key]}")
        if "metric" in chart and chart["metric"] not in metrics:
            fail(f"display.charts: unknown metric {chart['metric']}")
    if "map" in display:
        shown = display["map"].get("metric")
        if shown is not None and shown not in metrics:
            fail(f"display.map: unknown metric {shown}")
        for dim, member in display["map"].get("dims", {}).items():
            if dim not in allowed:
                fail(f"display.map: unknown dimension {dim}")
            elif member not in allowed[dim]:
                fail(f"display.map: {dim}={member} is not a member")
            elif shown in metrics and dim not in _metric_dims(metrics, shown):
                fail(f"display.map: {shown} is not keyed by {dim}")

    sets = vocab["boundary_sets"]
    sources = _sources(d)
    named = [s.get("boundary_set") for s in sources]
    if len(named) > 1 and (None in named or len(set(named)) != len(named)):
        fail("a list of sources must name a different boundary_set in each")
    for source in sources:
        if source["format"] == "long" and "metric" not in source:
            stored = [m for m, s in metrics.items() if "derive" not in s]
            templated = any("metric" in rule for rule in source.get("rules", []))
            if len(stored) != 1 and not templated:
                fail("a long source without a metric column needs rules that name the metric, or one stored metric")
        for column, dim in source.get("columns", {}).items():
            if dim not in used:
                fail(f"source.columns: {column} maps to unknown dimension {dim}")
        for rule in source.get("rules", []):
            for dim in rule.get("dims", {}):
                if dim not in used:
                    fail(f"source.rules: unknown dimension {dim}")
            if "metric" in rule and "{" not in rule["metric"] and rule["metric"] not in metrics:
                fail(f"source.rules: unknown metric {rule['metric']}")
        for column in source.get("value_columns", []):
            if "{" not in column["metric"] and column["metric"] not in metrics:
                fail(f"source.value_columns: unknown metric {column['metric']}")
            for dim in column.get("dims", {}):
                if dim not in used:
                    fail(f"source.value_columns: unknown dimension {dim}")
        boundary_set = source.get("boundary_set")
        if boundary_set is None:
            continue
        if boundary_set not in sets:
            fail(f"unknown boundary set {boundary_set}")
            continue
        kind = sets[boundary_set].get("kind", "admin")
        if kind != "admin" and "iso3" in source.get("keys", {}):
            fail(f"{boundary_set}: keys.iso3 is not used for {kind} sets")
        fixed = source.get("keys", {}).get("level")
        if isinstance(fixed, int) and fixed not in _level_numbers(sets[boundary_set]):
            fail(f"{boundary_set} has no level {fixed}")


def _metric_dims(metrics: dict, metric: str) -> list[str]:
    spec = metrics[metric]
    derive = spec.get("derive", {})
    base = derive.get("share_of") or derive.get("aggregate_of")
    return metrics[base].get("dimensions", []) if base in metrics else spec.get("dimensions", [])


def _build(vocab: dict, raw: dict[str, tuple[Path, dict]]) -> Catalogue:
    def opt(d: dict, key: str):
        return d.get(key)

    sections = {k: Section(k, v["title"], opt(v, "description")) for k, v in vocab["sections"].items()}
    units = {k: Unit(k, v["label"], v["quantity"]) for k, v in vocab["units"].items()}
    dimensions = {
        k: Dimension(
            k, v["label"],
            {mid: Member(mid, m["label"], opt(m, "short_label"), opt(m, "description"), opt(m, "value"), opt(m, "code"))
             for mid, m in _members_of(v).items()},
            opt(v, "description"), bool(v.get("ordered", False)), opt(v, "role"))
        for k, v in vocab["dimensions"].items()}
    boundary_sets = {
        k: BoundarySet(k, v["title"], v.get("kind", "admin"), tuple(v["levels"]), opt(v, "countries"),
                       opt(v, "attribution"), opt(v, "attribution_url"))
        for k, v in vocab["boundary_sets"].items()}
    section_order = list(sections)

    def sort_key(item):
        dataset_id, (_, d) = item
        return section_order.index(d["id"].split("/")[0]), d.get("order", 1 << 30), dataset_id

    datasets = {}
    for dataset_id, (path, d) in sorted(raw.items(), key=sort_key):
        metrics = {
            k: Metric(k, m["label"], m["unit"], tuple(m.get("dimensions", [])), m.get("aggregation", "none"),
                      opt(m, "weight"), opt(m, "description"), opt(m, "source_note"), opt(m, "derive"))
            for k, m in d["metrics"].items()}
        dims = {
            k: DatasetDimension(k, tuple(v.get("members") or dimensions[k].members), opt(v, "default"),
                                dict(v.get("labels", {})))
            for k, v in d.get("dimensions", {}).items()}
        datasets[dataset_id] = Dataset(
            dataset_id, d["title"], metrics, dims, tuple(_sources(d)), dict(d.get("display", {})),
            d.get("kind", "values"), opt(d, "x"), opt(d, "description"), opt(d, "documentation"),
            opt(d, "source_note"), opt(d, "order"), path)
    return Catalogue(sections, units, dimensions, boundary_sets, dict(vocab["countries"]), datasets)


# --- reading source files ----------------------------------------------------

@dataclass(frozen=True)
class Placed:
    """One value, placed in the cube."""
    metric: str
    dims: tuple[tuple[str, str], ...]             # sorted (dimension, member) pairs
    boundary_set: str | None                      # None: the country's admin set
    iso3: str | None
    level: int | None
    area_id: str
    value: float | None                           # None: published empty
    x: float | None = None                        # curve datasets

    @property
    def series(self) -> tuple[str, tuple[tuple[str, str], ...]]:
        return self.metric, self.dims


@dataclass(frozen=True)
class Skipped:
    reason: str                                   # the rule's `skip` value


@dataclass(frozen=True)
class Unplaced:
    """A row (or a file, when line is None) the rules cannot place. A loader must fail on these."""
    message: str


@dataclass(frozen=True)
class SourceRow:
    path: Path
    line: int | None                              # 1 is the header
    result: Placed | Skipped | Unplaced


def source_files(source: dict, data_dir: str | Path) -> list[tuple[Path, dict[str, str]]]:
    """The files a source's pattern matches, with the `{iso3}` and `{level}` they were matched with."""
    root = Path(data_dir)
    pattern = source["files"]
    regex, glob_pattern, pos = "", "", 0
    for found in re.finditer(r"\{(iso3|level)\}|\*", pattern):
        literal = pattern[pos:found.start()]
        regex += re.escape(literal)
        glob_pattern += literal
        name = found.group(1)
        if name is None:
            regex += "[^/]*"
        elif f"(?P<{name}>" in regex:
            regex += f"(?P={name})"
        else:
            regex += f"(?P<{name}>[^/]+)"
        glob_pattern += "*"
        pos = found.end()
    regex += re.escape(pattern[pos:])
    glob_pattern += pattern[pos:]
    out = []
    for path in sorted(root.glob(glob_pattern)):
        found = re.fullmatch(regex, path.relative_to(root).as_posix())
        if found:
            out.append((path, {k: v for k, v in found.groupdict().items() if v is not None}))
    return out


class _Rows:
    """Placement of rows of one source of one dataset."""

    def __init__(self, catalogue: Catalogue, dataset: Dataset, source: dict):
        self.catalogue, self.dataset, self.source = catalogue, dataset, source
        self.column_of = {dim: dim for dim in dataset.dimensions}
        self.column_of.update({dim: col for col, dim in source.get("columns", {}).items()})
        self.stored = list(dataset.stored_metrics)

    def fill(self, template: str, variables: dict, dim: str | None) -> str:
        def sub(match):
            value = variables[match.group(1)]
            if match.group(2) == "code":
                return self.catalogue.dimensions[dim].by_code(value)
            return value
        return _TEMPLATE.sub(sub, template)

    def member(self, dim: str, value: str) -> str:
        if value not in self.dataset.dimensions[dim].members:
            raise ValueError(f"{dim}={value!r} is not a member")
        return value

    def place_long(self, row: dict) -> tuple[str, dict] | str:
        """(metric, dims), or the skip reason. Raises ValueError/KeyError when the row cannot be placed."""
        source = self.source
        variables, metric, dims = dict(row), None, {}
        for rule in source.get("rules", []):
            when = rule.get("when", {})
            if not all(row.get(c) in ([v] if isinstance(v, str) else v) for c, v in when.items()):
                continue
            found = {c: re.fullmatch(p, row.get(c) or "") for c, p in rule.get("match", {}).items()}
            if not all(found.values()):
                continue
            if "skip" in rule:
                return rule["skip"]
            for m in found.values():
                variables.update({k: v for k, v in m.groupdict().items() if v is not None})
            metric = rule.get("metric")
            dims = dict(rule.get("dims", {}))
            break
        if metric is None:
            if "metric" in source:
                metric = row[source["metric"]]
            elif len(self.stored) == 1:
                metric = self.stored[0]
        else:
            metric = self.fill(metric, variables, None)
        if metric not in self.stored:
            raise ValueError(f"no stored metric for row (got {metric!r})")
        out = {}
        for dim in self.dataset.metrics[metric].dimensions:
            value = self.fill(dims[dim], variables, dim) if dim in dims else row.get(self.column_of[dim]) or ""
            if value == "":
                continue                          # the dimension does not apply to this series
            out[dim] = self.member(dim, value)
        unexpected = set(dims) - set(self.dataset.metrics[metric].dimensions)
        if unexpected:
            raise ValueError(f"{metric}: rule assigns {sorted(unexpected)}, which the metric is not keyed by")
        if "unit" in source and row[source["unit"]] != self.dataset.metrics[metric].unit:
            raise ValueError(f"{metric}: unit {row[source['unit']]!r}, expected {self.dataset.metrics[metric].unit!r}")
        return metric, out

    def known_columns(self) -> set[str]:
        source, keys = self.source, self.source.get("keys", {})
        known = {keys.get("area"), keys.get("iso3"), source.get("value", "value"), source.get("unit"),
                 source.get("metric"), source.get("x"), *source.get("ignore", []), *source.get("columns", {}),
                 *self.dataset.dimensions,
                 *(c for r in source.get("rules", []) for c in [*r.get("when", {}), *r.get("match", {})])}
        if isinstance(keys.get("level"), dict):
            known.add(keys["level"]["column"])
        return known - {None}

    def wide_columns(self, header: list[str]) -> tuple[dict[str, tuple[str, dict]], list[str], list[str]]:
        """(value column -> (metric, dims), columns ignored, columns matching several value_columns)."""
        source, keys = self.source, self.source.get("keys", {})
        structural = {source.get("x"), keys.get("area"), keys.get("iso3"), *self.dataset.dimensions}
        if isinstance(keys.get("level"), dict):
            structural.add(keys["level"]["column"])
        mapped, ignored, ambiguous = {}, [], []
        for column in header:
            if column in structural:
                continue
            hits = []
            for spec in source["value_columns"]:
                match = re.fullmatch(spec["match"], column) if "match" in spec else column == spec.get("column")
                if match:
                    hits.append((spec, match.groupdict() if hasattr(match, "groupdict") else {}))
            if not hits:
                ignored.append(column)
            elif len(hits) > 1:
                ambiguous.append(column)
            else:
                spec, variables = hits[0]
                metric = self.fill(spec["metric"], variables, None)
                dims = {dim: self.member(dim, self.fill(t, variables, dim)) for dim, t in spec.get("dims", {}).items()}
                mapped[column] = (metric, dims)
        return mapped, ignored, ambiguous


def _number(text: str | None) -> float | None:
    return None if text is None or text.strip() == "" else float(text)


def read_source(catalogue: Catalogue, dataset: Dataset | str, source: dict | int, data_dir: str | Path
                ) -> Iterator[SourceRow]:
    """Apply one source's rules to its files under `data_dir`. One SourceRow per value (long: per row)."""
    if isinstance(dataset, str):
        dataset = catalogue.datasets[dataset]
    if isinstance(source, int):
        source = dataset.sources[source]
    rows = _Rows(catalogue, dataset, source)
    keys = source.get("keys", {})
    boundary_set = source.get("boundary_set")
    for path, from_path in source_files(source, data_dir):
        with open(path, newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            header = list(reader.fieldnames or [])

            if source["format"] == "wide":
                try:
                    mapped, _, ambiguous = rows.wide_columns(header)
                except (KeyError, ValueError) as error:
                    yield SourceRow(path, 1, Unplaced(f"value columns: {error}"))
                    continue
                if ambiguous:
                    yield SourceRow(path, 1, Unplaced(f"columns match several value_columns: {ambiguous}"))
                    continue
                for metric, _ in mapped.values():
                    if metric not in rows.stored:
                        yield SourceRow(path, 1, Unplaced(f"no stored metric {metric!r}"))
                        break
                else:
                    row_dims = [d for d in dataset.dimensions if d in header]
                    for number, row in enumerate(reader, 2):
                        yield from _wide_row(rows, path, number, row, mapped, row_dims, keys, from_path, boundary_set)
                continue

            unknown = [c for c in header if c not in rows.known_columns()]
            if unknown:
                yield SourceRow(path, 1, Unplaced(f"unmapped columns {unknown}"))
                continue
            for number, row in enumerate(reader, 2):
                try:
                    placed = rows.place_long(row)
                    if isinstance(placed, str):
                        yield SourceRow(path, number, Skipped(placed))
                        continue
                    metric, dims = placed
                    where = _where(row, keys, from_path, boundary_set)
                    value = _number(row[source.get("value", "value")])
                except (KeyError, ValueError) as error:
                    yield SourceRow(path, number, Unplaced(f"{type(error).__name__}: {error}"))
                    continue
                yield SourceRow(path, number, Placed(metric, tuple(sorted(dims.items())), boundary_set, *where, value))


def _wide_row(rows: _Rows, path, number, row, mapped, row_dims, keys, from_path, boundary_set):
    try:
        where = _where(row, keys, from_path, boundary_set)
        extra = {dim: rows.member(dim, row[dim]) for dim in row_dims if row[dim] != ""}
        x = _number(row[rows.source["x"]]) if "x" in rows.source else None
    except (KeyError, ValueError) as error:
        yield SourceRow(path, number, Unplaced(f"{type(error).__name__}: {error}"))
        return
    for column, (metric, dims) in mapped.items():
        keyed = {**extra, **dims}
        try:
            unexpected = set(keyed) - set(rows.dataset.metrics[metric].dimensions)
            if unexpected:
                raise ValueError(f"{column}: {metric} is not keyed by {sorted(unexpected)}")
            value = _number(row[column])
        except (KeyError, ValueError) as error:
            yield SourceRow(path, number, Unplaced(f"{type(error).__name__}: {error}"))
            continue
        yield SourceRow(path, number, Placed(metric, tuple(sorted(keyed.items())), boundary_set, *where, value, x))


def _where(row: dict, keys: dict, from_path: dict, boundary_set: str | None) -> tuple[str | None, int | None, str]:
    """(iso3, level, area_id) of a row."""
    iso3 = row[keys["iso3"]] if "iso3" in keys else from_path.get("iso3")
    if "iso3" in keys and "iso3" in from_path and iso3 != from_path["iso3"]:
        raise ValueError(f"{keys['iso3']}={iso3!r} but the path is for {from_path['iso3']}")
    level = keys.get("level")
    if isinstance(level, dict):
        value = row[level["column"]]
        if value not in level["values"]:
            raise ValueError(f"{level['column']}={value!r} is not a known level")
        level = level["values"][value]
    elif level is None and "level" in from_path:
        level = int(from_path["level"])
    if "area" in keys:
        area = row[keys["area"]]
        if area == "":
            raise ValueError(f"empty area id in {keys['area']}")
    else:                                         # country-level values only
        if iso3 is None:
            raise ValueError("no area column and no iso3")
        area, level = iso3, 0 if level is None else level
    return iso3, level, area
