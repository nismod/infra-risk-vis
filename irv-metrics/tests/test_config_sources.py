"""The `source` rules of every GRI dataset place every row of every published CSV.

Skipped when the CSVs are absent. They are looked for in `csv-data/` next to
`backend/` (a symlink in a worktree), or in $CSV_DATA.

Run with `-s` to see the per-dataset counts.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
import os
from pathlib import Path

import pytest

from irv_metrics.config import CONFIG_DIR, Placed, Skipped, Unplaced, load_catalogue, read_source, source_files, _Rows

DATA = Path(os.environ.get("CSV_DATA") or Path(__file__).resolve().parents[2] / "csv-data")
pytestmark = pytest.mark.skipif(not DATA.is_dir(), reason=f"published CSVs not found at {DATA}")

# Stored metrics that the published CSVs do not carry; they come from the upstream API only.
NOT_IN_CSVS = {
    ("hazard/tropical-cyclone-wind", "area_in_band"),        # profile `bands` only; = differences of thresholds
    ("hazard/tropical-cyclone-wind", "area_in_band_share"),
}


@dataclass
class Report:
    files: int = 0
    rows: int = 0                                 # values: rows for long files, cells for wide ones
    placed: int = 0
    skipped: Counter = field(default_factory=Counter)
    empty: int = 0                                # placed with an empty value
    series: set = field(default_factory=set)
    unplaced: list = field(default_factory=list)
    duplicates: int = 0
    paths: set = field(default_factory=set)


@pytest.fixture(scope="module")
def catalogue():
    return load_catalogue()


@pytest.fixture(scope="module")
def reports(catalogue) -> dict[str, Report]:
    out = {}
    for dataset in catalogue.datasets.values():
        report = out[dataset.id] = Report()
        for source in dataset.sources:
            files = source_files(source, DATA)
            report.files += len(files)
            report.paths |= {path for path, _ in files}
            seen = Counter()
            for row in read_source(catalogue, dataset, source, DATA):
                result = row.result
                report.rows += 1
                if isinstance(result, Placed):
                    report.placed += 1
                    report.series.add(result.series)
                    report.empty += result.value is None
                    seen[(result.series, result.boundary_set, result.area_id, result.x)] += 1
                elif isinstance(result, Skipped):
                    report.skipped[result.reason] += 1
                else:
                    assert isinstance(result, Unplaced)
                    report.unplaced.append(f"{row.path.relative_to(DATA)}:{row.line}: {result.message}")
            report.duplicates += sum(n - 1 for n in seen.values() if n > 1)
    return out


def test_every_published_csv_is_read_by_one_source(catalogue, reports):
    claimed = Counter(path for report in reports.values() for path in report.paths)
    published = sorted(p for p in DATA.glob("[A-Z][A-Z][A-Z]/**/*.csv"))
    unread = [str(p.relative_to(DATA)) for p in published if p not in claimed]
    twice = [str(p.relative_to(DATA)) for p, n in claimed.items() if n > 1]
    assert not unread, f"published CSVs that no dataset reads: {unread}"
    assert not twice, f"CSVs read by several sources: {twice}"


DATASET_IDS = [p.relative_to(CONFIG_DIR).with_suffix("").as_posix()
               for p in sorted(CONFIG_DIR.rglob("*.yaml")) if not p.name.startswith("vocabulary")]


@pytest.mark.parametrize("dataset_id", DATASET_IDS)
def test_rules_place_every_row(catalogue, reports, dataset_id):
    report = reports[dataset_id]
    assert report.files, f"{dataset_id}: no source files"
    assert not report.unplaced, f"{dataset_id}: {len(report.unplaced)} unplaced, e.g. " + "\n".join(report.unplaced[:10])
    assert report.placed + sum(report.skipped.values()) == report.rows
    assert report.duplicates == 0, f"{dataset_id}: {report.duplicates} values placed on an already-filled (series, area)"
    dataset = catalogue.datasets[dataset_id]
    loaded = {metric for metric, _ in report.series}
    missing = set(dataset.stored_metrics) - loaded - {m for d, m in NOT_IN_CSVS if d == dataset_id}
    assert not missing, f"{dataset_id}: stored metrics with no values in the CSVs: {missing}"


def test_wide_sources_map_every_column(catalogue):
    import csv
    for dataset in catalogue.datasets.values():
        for source in dataset.sources:
            if source["format"] != "wide":
                continue
            for path, _ in source_files(source, DATA):
                with open(path, newline="") as handle:
                    header = next(csv.reader(handle))
                mapped, ignored, ambiguous = _Rows(catalogue, dataset, source).wide_columns(header)
                assert not ignored and not ambiguous, (path, ignored, ambiguous)


def test_report(reports, capsys):
    lines = [f"{'dataset':52} {'files':>5} {'rows':>8} {'placed':>8} {'skipped':>8} {'empty':>6} {'series':>6}"]
    for dataset_id, r in reports.items():
        skipped = sum(r.skipped.values())
        reasons = ", ".join(f"{k} {v}" for k, v in r.skipped.items())
        lines.append(f"{dataset_id:52} {r.files:>5} {r.rows:>8} {r.placed:>8} {skipped:>8} {r.empty:>6} {len(r.series):>6}"
                     + (f"  ({reasons})" if reasons else ""))
    total = [sum(getattr(r, k) for r in reports.values()) for k in ("files", "rows", "placed")]
    lines.append(f"{'total':52} {total[0]:>5} {total[1]:>8} {total[2]:>8} "
                 f"{sum(sum(r.skipped.values()) for r in reports.values()):>8} "
                 f"{sum(r.empty for r in reports.values()):>6} {sum(len(r.series) for r in reports.values()):>6}")
    with capsys.disabled():
        print("\n" + "\n".join(lines))
