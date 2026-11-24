"""Load published metric CSVs and boundary files into the cube, driven by the YAML dataset config.

    python -m irv_metrics.loaders.csv --data-dir DIR [--config backend/config/gri] [--countries KEN MOZ]
                                      [--if-loaded skip|replace] [--reset] [--report out.json]

File patterns in the config (dataset `source.files`, boundary set `levels[].files`)
are relative to the data directory. For GRI on the dev machine:

    data/
      KEN -> .../csv-data/KEN                                published CSVs, one directory per country
      MOZ -> .../csv-data/MOZ
      boundaries/cgaz -> /data/projects/svi/incoming_data/boundaries
                                                             geoBoundariesCGAZ_ADM{0,1,2}.gpkg (global files)

Discovery. Countries and levels come from the files each source's pattern
matches (`{iso3}`, `{level}`). As in the original ETL:
- an admin-0 file without a single non-blank value is a blank stub and is not
  loaded (so it does not become all-null country rows);
- a (dataset, country) with admin-1/2 files but no usable admin-0 file fails the
  load: country totals are never made up from subnational rows;
- at admin level 0 an area id that is the country's name instead of its ISO3
  (some drops write `Kenya`) is read as the ISO3.

Order, per country, in one transaction: the country row; its areas in its admin
set, read from the set's level files (only that country's features); parents
from the `parent` property or by location; overlay sets that take their
countries from the set are re-related; then each dataset: remove its old
values, a load record per (dataset, boundary set, country), COPY the values (and
curve points). Overlay and group sets are loaded globally (iso3 NULL) after the
countries; their areas only when the set has none yet (`--reset` reloads them).

Structure. The cube's sections, units, dimensions, boundary sets, datasets
(under their config ids), metrics and members are seeded from the config
(`common.seed_catalogue`) before anything is loaded.

Values computed at load. Some stored metrics are not published but follow from
published ones; `COMPUTED` declares them per dataset id (see `Bands`). For GRI:
the tropical-cyclone wind bands `area_in_band[_share]` are the differences of
consecutive `area_exceeding[_share]` thresholds.

Countries. A country is loaded when a file pattern with `{iso3}` matches its
files, or, when some admin-set source reads files without `{iso3}` (one file
for a single-country deployment), when the vocabulary lists it under
`countries`.

A row the source rules cannot place, a value on an unknown area or a value
placed twice fails the country, which is rolled back; the loader goes on with
the next and exits 1. Published rows without a value are stored as NULL.

DATABASE_URL selects the database (default postgresql://gri:gri@localhost:5433/gri).
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass, field
import json
from pathlib import Path
import sys
import tempfile
import time

from ..config import CONFIG_DIR, Catalogue, Dataset, Placed, Skipped, Unplaced, load_catalogue, read_source, source_files
from .common import (SEARCH_PATH, SchemaManagedError, apply_schema, begin_load, check_series, connect, country_loaded,
                     ensure_series, finish_load, insert_area, log, relate_overlay, remove_country, remove_observations,
                     seed_catalogue, set_parents_by_location)


class LoadError(Exception):
    def __init__(self, problems: list[str]):
        self.problems = problems
        super().__init__(f"{len(problems)} problem(s): " + "; ".join(problems[:5]))


# --- discovery -----------------------------------------------------------------

@dataclass(frozen=True)
class SourceFile:
    dataset: str
    source: int                                   # index into the dataset's sources
    path: Path
    rel: str                                      # relative to the data dir
    iso3: str | None                              # from the path; None: the file may hold several countries
    level: int | None                             # from the path


@dataclass
class Discovery:
    files: list[SourceFile] = field(default_factory=list)     # usable
    stubs: list[SourceFile] = field(default_factory=list)     # blank admin-0 files, not loaded
    gaps: list[str] = field(default_factory=list)             # admin-0 completeness failures
    configured: list[str] = field(default_factory=list)       # vocabulary countries, for files without {iso3}

    def countries(self) -> list[str]:
        return sorted({f.iso3 for f in self.files if f.iso3} | set(self.configured))

    def of(self, iso3: str | None) -> list[SourceFile]:
        return [f for f in self.files if f.iso3 == iso3]


def has_values(path: Path, column: str) -> bool:
    """True when at least one row has a non-blank value."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return any((row.get(column) or "").strip() for row in csv.DictReader(handle))


def source_set(catalogue: Catalogue, source: dict) -> str | None:
    """The source's boundary set; None for the country's admin set."""
    return source.get("boundary_set")


def is_admin(catalogue: Catalogue, set_id: str | None) -> bool:
    return set_id is None or catalogue.boundary_sets[set_id].kind == "admin"


def discover(catalogue: Catalogue, data_dir: Path) -> Discovery:
    found = Discovery()
    for dataset in catalogue.datasets.values():
        for index, source in enumerate(dataset.sources):
            levels: dict[str, set[int]] = defaultdict(set)
            for path, variables in source_files(source, data_dir):
                level = int(variables["level"]) if "level" in variables else None
                entry = SourceFile(dataset.id, index, path, path.relative_to(data_dir).as_posix(),
                                   variables.get("iso3"), level)
                if (level == 0 and source["format"] == "long"
                        and not has_values(path, source.get("value", "value"))):
                    found.stubs.append(entry)
                    continue
                found.files.append(entry)
                if entry.iso3 and level is not None:
                    levels[entry.iso3].add(level)
            if not is_admin(catalogue, source_set(catalogue, source)):
                continue
            for iso3, present in sorted(levels.items()):
                if any(level > 0 for level in present) and 0 not in present:
                    found.gaps.append(f"{iso3} {dataset.id}: admin levels {sorted(present)} without usable admin 0")
    if any(f.iso3 is None and is_admin(catalogue, source_set(catalogue, catalogue.datasets[f.dataset].sources[f.source]))
           for f in found.files):
        found.configured = list(catalogue.countries)
    return found


def admin_set_for(catalogue: Catalogue, iso3: str, files: list[SourceFile]) -> str:
    """The vocabulary's `countries.<iso3>.boundary_set`; else the first admin set a source of the
    country's files names; else the first admin set in the vocabulary."""
    configured = catalogue.countries.get(iso3, {}).get("boundary_set")
    if configured:
        return configured
    for f in files:
        named = catalogue.datasets[f.dataset].sources[f.source].get("boundary_set")
        if named and catalogue.boundary_sets[named].kind == "admin":
            return named
    for set_id, spec in catalogue.boundary_sets.items():
        if spec.kind == "admin":
            return set_id
    raise LoadError([f"{iso3}: the vocabulary has no admin boundary set"])


# --- boundaries ----------------------------------------------------------------

@dataclass(frozen=True)
class Area:
    level: int
    area_id: str
    name: str
    geometry: bytes | None                        # WKB
    parent_id: str | None
    iso3: str | None                              # admin: the area's country; group: a member country
    source: str


def read_level(spec_level: dict, number: int, kind: str, data_dir: Path, iso3: str | None) -> list[Area]:
    """The features of one level file; for an admin set, only those of `iso3`.

    GeoPackages and other OGR formats are read with pyogrio; `.csv` files (group
    membership tables) with the csv module. Where the layer has a geoBoundaries
    `shapeType` field, only features of this level (`ADM{n}`) are kept: the CGAZ
    ADM2 file also carries some coarser polygons.
    """
    pattern = spec_level["files"]
    if "{iso3}" in pattern and iso3 is None:
        raise LoadError([f"{pattern}: {{iso3}} in a {kind} set's level file"])
    rel = pattern.replace("{iso3}", iso3 or "")
    path = data_dir / rel
    if not path.exists():
        raise LoadError([f"boundary file {path} not found"])
    id_col, name_col = spec_level["id"], spec_level["name"]
    iso3_col, parent_col = spec_level.get("iso3"), spec_level.get("parent")
    by_column = kind == "admin" and iso3 is not None and iso3_col is not None and "{iso3}" not in pattern

    if path.suffix.lower() == ".csv":
        with open(path, newline="", encoding="utf-8-sig") as handle:
            records = list(csv.DictReader(handle))
        if by_column:
            records = [r for r in records if r.get(iso3_col) == iso3]
        geometries = [None] * len(records)
    else:
        import pyogrio
        import pyogrio.raw

        fields = list(pyogrio.read_info(path)["fields"])
        columns = [c for c in dict.fromkeys([id_col, name_col, iso3_col, parent_col]) if c]
        missing = [c for c in columns if c not in fields]
        if missing:
            raise LoadError([f"{rel}: no field(s) {missing} (has {fields})"])
        if "shapeType" in fields:
            columns.append("shapeType")
        where = None
        if by_column:
            if not iso3.isalnum():
                raise LoadError([f"bad ISO3 {iso3!r}"])
            where = f"\"{iso3_col}\" = '{iso3}'"
        meta, _, geometries, data = pyogrio.raw.read(path, columns=columns, where=where)
        names = list(meta["fields"])
        records = [dict(zip(names, values)) for values in zip(*data)] if len(names) else [{} for _ in geometries]
        if "shapeType" in names:
            keep = [i for i, r in enumerate(records) if r["shapeType"] == f"ADM{number}"]
            records, geometries = [records[i] for i in keep], [geometries[i] for i in keep]

    out = []
    for record, geometry in zip(records, geometries):
        area_id = str(record[id_col] or "")
        if not area_id:
            raise LoadError([f"{rel}: a feature without {id_col}"])
        if kind != "group" and geometry is None:
            raise LoadError([f"{rel}: {area_id} has no geometry"])
        name = record.get(name_col) or area_id      # a few published shapes have no name
        member = record.get(iso3_col) if iso3_col else iso3
        out.append(Area(number, area_id, str(name), geometry, record.get(parent_col) if parent_col else None,
                        member, rel))
    return out


def read_boundaries(catalogue: Catalogue, set_id: str, data_dir: Path, iso3: str | None = None) -> list[Area]:
    spec = catalogue.boundary_sets[set_id]
    areas = []
    for position, level in enumerate(spec.levels):
        areas += read_level(level, level.get("level", position), spec.kind, data_dir, iso3)
    return areas


def insert_areas(cur, set_id: str, areas: list[Area]) -> dict[str, int]:
    keys = {}
    for area in sorted(areas, key=lambda a: a.level):
        if area.area_id in keys:
            if area.geometry is None:               # group tables: one row per (group, member)
                continue
            raise LoadError([f"{set_id}: area {area.area_id} appears twice"])
        keys[area.area_id] = insert_area(cur, set_id, area.area_id, area.level, area.name, area.geometry,
                                         source=area.source)
    return keys


def set_parents(cur, set_id: str, areas: list[Area], keys: dict[str, int], iso3: str | None):
    """Parents from the level files' `parent` property; levels without one by location."""
    levels = sorted({a.level for a in areas})
    by_location = []
    for level in levels[1:]:
        at_level = [a for a in areas if a.level == level]
        if all(a.parent_id is None for a in at_level):
            by_location.append(level)
            continue
        unknown = [a.parent_id for a in at_level if a.parent_id not in keys]
        if unknown:
            raise LoadError([f"{set_id} level {level}: unknown parent ids, e.g. {unknown[:3]}"])
        cur.executemany("UPDATE area SET parent_area_key = %s WHERE area_key = %s",
                        [(keys[a.parent_id], keys[a.area_id]) for a in at_level])
    for level in by_location:
        set_parents_by_location(cur, set_id, level, iso3)


def ensure_country(cur, iso3: str, name: str, admin_set: str | None):
    cur.execute("INSERT INTO country VALUES (%s, %s, %s) ON CONFLICT (iso3) DO UPDATE SET name = excluded.name,"
                " boundary_set = coalesce(country.boundary_set, excluded.boundary_set)", (iso3, name, admin_set))


def area_index(cur, set_id: str, iso3: str | None) -> dict[str, tuple[int, int]]:
    """area_id -> (area_key, level) for a set's areas (related to iso3, if given)."""
    cur.execute("SELECT area_id, area_key, level FROM area WHERE boundary_set = %s"
                " AND (%s::text IS NULL OR area_key IN (SELECT area_key FROM area_country WHERE iso3 = %s))",
                (set_id, iso3, iso3))
    return {area_id: (key, level) for area_id, key, level in cur.fetchall()}


# --- values ----------------------------------------------------------------------

@dataclass
class DatasetStats:
    files: int = 0
    rows: int = 0                                 # values read: rows of long files, cells of wide ones
    placed: int = 0
    skipped: Counter = field(default_factory=Counter)
    empty: int = 0                                # placed without a value: stored as NULL (curves: dropped)
    other_country: int = 0                        # rows of another country in a multi-country file
    computed: int = 0                             # values computed at load (COMPUTED)
    series: int = 0
    seconds: float = 0.0

    def line(self) -> str:
        skipped = sum(self.skipped.values())
        reasons = f" ({', '.join(f'{k} {v}' for k, v in sorted(self.skipped.items()))})" if skipped else ""
        other = f", other countries {self.other_country}" if self.other_country else ""
        computed = f", computed {self.computed}" if self.computed else ""
        return (f"files {self.files}, rows {self.rows}, placed {self.placed}, skipped {skipped}{reasons},"
                f" empty {self.empty}{other}{computed}, series {self.series}, {self.seconds:.1f}s")


# --- values computed at load ----------------------------------------------------

Cells = dict[tuple, dict[int, float | None]]      # (metric, dims) -> area_key -> value


@dataclass(frozen=True)
class Bands:
    """Band values from cumulative thresholds: band i = threshold i - threshold i+1; the top band is the
    top threshold.

    `metrics` maps each band metric to the threshold metric it is computed from.
    Threshold and band members pair up in the dataset's member order (the first
    threshold, "at or above the lowest class", with the lowest band). A band
    value is NULL when either threshold is NULL, and absent when either is missing.
    The band metrics are stored metrics without source rules; a source that
    publishes them too is an error.
    """
    metrics: dict[str, str]
    threshold: str                                # the thresholds' dimension
    band: str                                     # the bands' dimension

    def check(self, dataset: Dataset) -> list[str]:
        problems = []
        thresholds = dataset.dimensions.get(self.threshold)
        bands = dataset.dimensions.get(self.band)
        if thresholds is None or bands is None:
            return [f"{dataset.id}: needs dimensions {self.threshold} and {self.band}"]
        if len(thresholds.members) != len(bands.members):
            problems.append(f"{dataset.id}: {len(thresholds.members)} {self.threshold} members,"
                            f" {len(bands.members)} {self.band} members")
        for band_metric, threshold_metric in self.metrics.items():
            band, threshold = dataset.stored_metrics.get(band_metric), dataset.stored_metrics.get(threshold_metric)
            if band is None or threshold is None:
                problems.append(f"{dataset.id}: {band_metric} and {threshold_metric} must be stored metrics")
                continue
            expected = tuple(self.band if d == self.threshold else d for d in threshold.dimensions)
            if self.threshold not in threshold.dimensions or band.dimensions != expected:
                problems.append(f"{dataset.id}: {band_metric} must be keyed like {threshold_metric},"
                                f" with {self.band} for {self.threshold} ({list(expected)})")
            if band.unit != threshold.unit:
                problems.append(f"{dataset.id}: {band_metric} and {threshold_metric} differ in unit")
        return problems

    def compute(self, dataset: Dataset, values: Cells) -> Cells:
        thresholds = dataset.dimensions[self.threshold].members
        bands = dict(zip(thresholds, dataset.dimensions[self.band].members))
        upper = dict(zip(thresholds, thresholds[1:]))           # the next threshold up; none for the top one
        out: Cells = {}
        for band_metric, threshold_metric in self.metrics.items():
            for (metric, dims), cells in values.items():
                if metric != threshold_metric:
                    continue
                keyed = dict(dims)
                lower = keyed.get(self.threshold)
                if lower not in bands:
                    continue
                above = values.get((metric, tuple(sorted({**keyed, self.threshold: upper[lower]}.items()))), {}) \
                    if lower in upper else None
                band = {}
                for key, low in cells.items():
                    if above is not None and key not in above:
                        continue
                    high = 0.0 if above is None else above[key]
                    band[key] = None if low is None or high is None else low - high
                if band:
                    dims_out = {k: v for k, v in keyed.items() if k != self.threshold} | {self.band: bands[lower]}
                    out[(band_metric, tuple(sorted(dims_out.items())))] = band
        return out


# Stored metrics computed at load, by dataset id. C2 has no construct for these (see the README).
COMPUTED: dict[str, tuple[Bands, ...]] = {
    "hazard/tropical-cyclone-wind": (
        Bands({"area_in_band": "area_exceeding", "area_in_band_share": "area_exceeding_share"},
              threshold="wind_threshold", band="wind_band"),
    ),
}


def computed_problems(catalogue: Catalogue) -> list[str]:
    """Declarations in COMPUTED that the catalogue's datasets do not fit (datasets it lacks are ignored)."""
    return [p for dataset_id, rules in COMPUTED.items() if dataset_id in catalogue.datasets
            for rule in rules for p in rule.check(catalogue.datasets[dataset_id])]


def load_dataset(cur, catalogue: Catalogue, dataset_id: str, source_index: int, set_id: str,
                 iso3: str | None, read_dir: Path, aliases: set[str] = frozenset(),
                 data_dir: Path | None = None) -> DatasetStats:
    """Read one source of one dataset under `read_dir` and replace its values on (set_id, iso3)."""
    started = time.monotonic()
    dataset = catalogue.datasets[dataset_id]
    source = dataset.sources[source_index]
    curve = dataset.kind == "curve"
    areas = area_index(cur, set_id, iso3)
    stats = DatasetStats(files=len(source_files(source, read_dir)))
    problems: list[str] = []
    unknown: Counter = Counter()
    values: dict[tuple, dict[int, float | None]] = defaultdict(dict)      # (metric, dims) -> area_key -> value
    points: dict[tuple, dict[int, list]] = defaultdict(lambda: defaultdict(list))
    multi = "{iso3}" not in source["files"]

    for row in read_source(catalogue, dataset, source, read_dir):
        stats.rows += 1
        result = row.result
        where = f"{row.path.relative_to(read_dir).as_posix()}:{row.line}"
        if isinstance(result, Skipped):
            stats.skipped[result.reason] += 1
            continue
        if isinstance(result, Unplaced):
            problems.append(f"{where}: {result.message}")
            continue
        assert isinstance(result, Placed)
        if iso3 is not None and result.iso3 is not None and result.iso3 != iso3:
            if multi:
                stats.other_country += 1
                continue
            problems.append(f"{where}: a row for {result.iso3} in a file for {iso3}")
            continue
        area_id = result.area_id
        if area_id not in areas and result.level in (0, None) and area_id in aliases:
            area_id = iso3
        found = areas.get(area_id)
        if found is None:
            unknown[area_id] += 1
            continue
        key, level = found
        if result.level is not None and result.level != level:
            problems.append(f"{where}: {area_id} is at level {level}, the row says {result.level}")
            continue
        stats.placed += 1
        if curve:
            if result.value is None:
                stats.empty += 1
                continue
            points[result.series][key].append((result.x, result.value))
            continue
        stats.empty += result.value is None
        cell = values[result.series]
        if key in cell:
            problems.append(f"{where}: {result.metric} {dict(result.dims)} on {area_id} placed twice")
            continue
        cell[key] = result.value
    if unknown:
        problems.append(f"{sum(unknown.values())} values on {len(unknown)} areas not in {set_id}"
                        f"{' for ' + iso3 if iso3 else ''}, e.g. {list(unknown)[:5]}")
    if not curve:
        for rule in COMPUTED.get(dataset_id, ()):
            for key, cells in rule.compute(dataset, values).items():
                if key in values:
                    problems.append(f"{key[0]} {dict(key[1])} is published and computed at load")
                    continue
                values[key] = cells
                stats.computed += len(cells)
    if problems:
        raise LoadError([f"{dataset_id}: {p}" for p in problems])

    series_keys = set(values) | set(points)
    stats.series = len(series_keys)
    problems = check_series(catalogue, {(dataset_id, metric, dims) for metric, dims in series_keys})
    if problems:
        raise LoadError(problems)
    series = ensure_series(cur, {(dataset_id, metric, dims) for metric, dims in series_keys})
    remove_observations(cur, set_id, iso3, dataset_id)
    label = source["files"] if data_dir is None else (data_dir / source["files"]).as_posix()
    notes = {"files": stats.files, "rows": stats.rows, "placed": stats.placed, "skipped": dict(stats.skipped),
             "empty": stats.empty, "computed": stats.computed}
    load_id = begin_load(cur, dataset_id, set_id, iso3, f"csv {label}", notes)
    with cur.copy("COPY observation (series_id, area_key, value) FROM STDIN") as copy:
        for (metric, dims), cells in values.items():
            series_id = series[(dataset_id, metric, dims)]
            for key, value in cells.items():
                copy.write_row((series_id, key, value))
    with cur.copy("COPY curve_point (series_id, area_key, seq, x, y) FROM STDIN") as copy:
        for (metric, dims), by_area in points.items():
            series_id = series[(dataset_id, metric, dims)]
            for key, xy in by_area.items():
                for seq, (x, y) in enumerate(xy):
                    copy.write_row((series_id, key, seq, x, y))
    finish_load(cur, load_id)
    stats.seconds = time.monotonic() - started
    return stats


# --- the load ----------------------------------------------------------------------

@dataclass
class Loader:
    conn: object
    catalogue: Catalogue
    data_dir: Path
    if_loaded: str = "skip"
    commit: bool = True                           # False: leave the transaction to the caller (tests)
    report: dict = field(default_factory=dict)
    failed: list[str] = field(default_factory=list)

    def __post_init__(self):
        problems = computed_problems(self.catalogue)
        if problems:
            raise LoadError(problems)

    def prepare(self, reset: bool = False):
        with self.conn.cursor() as cur:
            apply_schema(cur, reset=reset)
            seed_catalogue(cur, self.catalogue)
        self._commit()

    def _commit(self):
        if self.commit:
            self.conn.commit()

    def run(self, countries: list[str] | None = None, reset: bool = False) -> bool:
        self.prepare(reset)
        found = discover(self.catalogue, self.data_dir)
        selected = countries or found.countries()
        for stub in found.stubs:
            if stub.iso3 in selected:
                log(f"skipped blank admin-0 stub: {stub.rel}")
        gaps = [g for g in found.gaps if g.split()[0] in selected]
        if gaps:
            for gap in gaps:
                log(f"error: {gap}")
            log("error: usable admin 1/2 without usable admin 0; publish the country rows or remove the stubs")
            return False
        missing = [iso3 for iso3 in selected if iso3 not in found.countries()]
        if missing:
            log(f"error: no files for {missing} under {self.data_dir}")
            return False
        for iso3 in selected:
            self.load_country(iso3, found)
        self.load_global_sets(found)
        if self.commit:
            with self.conn.cursor() as cur:
                cur.execute("ANALYZE")
            self.conn.commit()
        return not self.failed

    def load_country(self, iso3: str, found: Discovery):
        files = found.of(iso3) + [f for f in found.of(None)
                                  if is_admin(self.catalogue, self.catalogue.datasets[f.dataset].sources[f.source]
                                              .get("boundary_set"))]
        admin_set = admin_set_for(self.catalogue, iso3, files)
        started = time.monotonic()
        report = self.report[iso3] = {"boundary_set": admin_set, "datasets": {}}
        try:
            with self.conn.cursor() as cur:
                cur.execute(SEARCH_PATH)
                if country_loaded(cur, iso3, admin_set):
                    if self.if_loaded == "skip":
                        log(f"{iso3}: already loaded, skipping (use --if-loaded replace to reload)")
                        report["skipped"] = True
                        return
                    log(f"{iso3}: removing previous load")
                    remove_country(cur, iso3, admin_set)
                self._load_country(cur, iso3, admin_set, files, report)
            self._commit()
            log(f"{iso3}: loaded in {time.monotonic() - started:.1f}s")
        except LoadError as error:
            if self.commit:
                self.conn.rollback()
            self.failed.append(iso3)
            report["errors"] = error.problems
            log(f"{iso3}: FAILED, rolled back")
            for problem in error.problems[:20]:
                log(f"  {problem}")

    def _load_country(self, cur, iso3: str, admin_set: str, files: list[SourceFile], report: dict):
        started = time.monotonic()
        areas = read_boundaries(self.catalogue, admin_set, self.data_dir, iso3)
        level0 = [a for a in areas if a.level == 0]
        if [a.area_id for a in level0] != [iso3]:
            raise LoadError([f"{admin_set} level 0 for {iso3}: expected one area {iso3}, found {[a.area_id for a in level0]}"])
        configured = self.catalogue.countries.get(iso3, {}).get("name")
        name = configured or level0[0].name or iso3
        ensure_country(cur, iso3, name, admin_set)
        keys = insert_areas(cur, admin_set, areas)
        cur.executemany("INSERT INTO area_country VALUES (%s, %s, 1)", [(key, iso3) for key in keys.values()])
        set_parents(cur, admin_set, areas, keys, iso3)
        counts = Counter(a.level for a in areas)
        log(f"{iso3}: {admin_set}: " + ", ".join(f"level {k} {v}" for k, v in sorted(counts.items()))
            + f" areas in {time.monotonic() - started:.1f}s")
        report["areas"] = dict(counts)
        for set_id, spec in self.catalogue.boundary_sets.items():
            if spec.kind == "overlay" and (spec.countries or {}).get("from") == admin_set and area_index(cur, set_id, None):
                log(f"{iso3}: related to {relate_overlay(cur, set_id, iso3)} areas of {set_id}")
        aliases = {iso3, name, *(a.name for a in level0)}

        with tempfile.TemporaryDirectory(prefix=f"irv-csv-{iso3}-") as view:
            # Only this country's files, at their places under the data dir, so the source
            # patterns (and their {iso3}/{level} captures) read nothing else.
            view = Path(view)
            for f in files:
                (view / f.rel).parent.mkdir(parents=True, exist_ok=True)
                (view / f.rel).symlink_to(f.path.resolve())
            for dataset_id, index in sorted({(f.dataset, f.source) for f in files},
                                            key=lambda k: list(self.catalogue.datasets).index(k[0])):
                stats = load_dataset(cur, self.catalogue, dataset_id, index, admin_set, iso3, view, aliases,
                                     self.data_dir)
                report["datasets"][dataset_id] = stats.__dict__ | {"skipped": dict(stats.skipped)}
                log(f"{iso3}: {dataset_id:46} {stats.line()}")

    def load_global_sets(self, found: Discovery):
        """Overlay and group sets: areas when the set has none, then each dataset with sources on them."""
        for set_id, spec in self.catalogue.boundary_sets.items():
            if spec.kind == "admin":
                continue
            sources = [(d.id, i) for d in self.catalogue.datasets.values() for i, s in enumerate(d.sources)
                       if s.get("boundary_set") == set_id]
            if not sources:
                continue
            report = self.report[set_id] = {"datasets": {}}
            try:
                with self.conn.cursor() as cur:
                    cur.execute(SEARCH_PATH)
                    if not area_index(cur, set_id, None):
                        self._load_global_areas(cur, set_id, spec)
                    for dataset_id, index in sources:
                        cur.execute("SELECT 1 FROM load_record WHERE dataset_id = %s AND boundary_set = %s"
                                    " AND iso3 IS NULL", (dataset_id, set_id))
                        if cur.fetchone() and self.if_loaded == "skip":
                            log(f"{set_id}: {dataset_id} already loaded, skipping")
                            continue
                        stats = load_dataset(cur, self.catalogue, dataset_id, index, set_id, None, self.data_dir,
                                             data_dir=self.data_dir)
                        report["datasets"][dataset_id] = stats.__dict__ | {"skipped": dict(stats.skipped)}
                        log(f"{set_id}: {dataset_id:46} {stats.line()}")
                self._commit()
            except LoadError as error:
                if self.commit:
                    self.conn.rollback()
                self.failed.append(set_id)
                report["errors"] = error.problems
                log(f"{set_id}: FAILED, rolled back")
                for problem in error.problems[:20]:
                    log(f"  {problem}")

    def _load_global_areas(self, cur, set_id: str, spec):
        started = time.monotonic()
        areas = read_boundaries(self.catalogue, set_id, self.data_dir)
        keys = insert_areas(cur, set_id, areas)
        if spec.kind == "overlay":
            related = relate_overlay(cur, set_id)
            set_parents(cur, set_id, areas, keys, None)
        else:                                       # group: one row per (group, member country)
            cur.execute("SELECT iso3 FROM country")
            known = {row[0] for row in cur.fetchall()}
            pairs = {(keys[a.area_id], a.iso3) for a in areas if a.iso3 in known}
            cur.executemany("INSERT INTO area_country VALUES (%s, %s, 1)", sorted(pairs))
            related = len(pairs)
            unknown = sorted({a.iso3 for a in areas if a.iso3 not in known})
            if unknown:
                log(f"{set_id}: members not loaded as countries, not related: {unknown[:10]}")
        log(f"{set_id}: {len(keys)} areas, {related} country relations in {time.monotonic() - started:.1f}s")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path, required=True, help="root the config's file patterns are relative to")
    parser.add_argument("--config", type=Path, default=CONFIG_DIR, help=f"dataset config directory (default {CONFIG_DIR})")
    parser.add_argument("--countries", nargs="+", metavar="ISO3", help="load only these (default: every country found)")
    parser.add_argument("--if-loaded", choices=("skip", "replace"), default="skip",
                        help="what to do with a country (or global set dataset) that is already loaded (default: skip)")
    parser.add_argument("--reset", action="store_true", help="drop and recreate the metrics schema first (refused where Alembic manages the schema,"
                        " unless IRV_METRICS_ALLOW_RESET=1)")
    parser.add_argument("--report", type=Path, help="write a JSON load report")
    args = parser.parse_args(argv)
    if not args.data_dir.is_dir():
        parser.error(f"{args.data_dir} is not a directory")
    catalogue = load_catalogue(args.config)
    started = time.monotonic()
    with connect() as conn:
        loader = Loader(conn, catalogue, args.data_dir, args.if_loaded)
        try:
            ok = loader.run([c.upper() for c in args.countries] if args.countries else None, reset=args.reset)
        except SchemaManagedError as error:
            log(f"error: {error}")
            return 2
    log(f"{'done' if ok else 'FAILED: ' + ', '.join(loader.failed) if loader.failed else 'FAILED'}"
        f" in {time.monotonic() - started:.0f}s")
    if args.report:
        args.report.write_text(json.dumps(loader.report, indent=2, default=str) + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
