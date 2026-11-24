"""Map infra-risk-vis pipeline layers to irv-datapkg catalogue items.

A pipeline opts in to sourcing its raw rasters from the irv-datapkg STAC
catalogue with a ``pipelines/{pipeline}/datapkg.yml`` mapping file::

    dataset: aqueduct_flood        # catalogue dataset id (irv-datapkg docs/datasets.csv)
    properties:                    # layers.csv column -> STAC item property
      hazard: hazard
      rp: return_period
    values:                        # optional layers.csv value -> property value, per column
      hazard: {fluvial: river}
      gcm: {None: null}            # null: the property must be absent or null
    filters: {}                    # optional extra property filters for every row
    filename: "{id}.tif"           # new rows only: filename template
    columns: {var: "depth_{hazard}"}  # new rows only: other column templates

Each layers.csv row resolves to exactly one catalogue item of ``dataset`` (or
of the row's optional ``datapkg_dataset`` column, for pipelines such as
``nature`` that combine several catalogue datasets), by (in order):

1. ``datapkg_item`` column: the STAC item id;
2. otherwise the property filter built from ``properties``/``values``/``filters``,
   updated with the row's optional ``datapkg_filter`` column
   (``key=value;key=v1|v2``, a value ``null`` meaning absent).

This module is imported by the Snakemake rules (no irv_datapkg import at
module level, so parsing the workflow does not need the client installed) and
by the ``datapkg_fetch``, ``layers_from_catalogue`` and
``metadata_from_catalogue`` scripts.
"""

import json
import math
import os
import string
from pathlib import Path
from typing import Optional

import pandas as pd
import yaml

ITEM_COLUMN = "datapkg_item"
FILTER_COLUMN = "datapkg_filter"
DATASET_COLUMN = "datapkg_dataset"
NULL = None  # in filters: property absent or null (the client's Catalog.ABSENT)


class MappingError(Exception):
    pass


# -- configuration ----------------------------------------------------------


def datapkg_config(config: dict) -> dict:
    """The ``datapkg`` section of the workflow config, with defaults."""
    cfg = dict(
        base=None,
        private_base=None,
        release=None,
        lock="datapkg.lock.yaml",
        cache_dir=None,
        enabled=[],
    )
    cfg.update((config or {}).get("datapkg") or {})
    cfg["enabled"] = list(cfg.get("enabled") or [])
    return cfg


def pipeline_dir(pipeline: str, root: str = ".") -> Path:
    return Path(root, "pipelines", pipeline)


def load_mapping(pipeline: str, root: str = ".") -> dict:
    path = pipeline_dir(pipeline, root) / "datapkg.yml"
    if not path.exists():
        raise MappingError(
            f"Pipeline {pipeline} is enabled for irv-datapkg but has no {path}"
        )
    with open(path) as f:
        mapping = yaml.safe_load(f) or {}
    if not mapping.get("dataset"):
        raise MappingError(f"{path} must name the catalogue 'dataset'")
    mapping.setdefault("properties", {})
    mapping.setdefault("values", {})
    mapping.setdefault("filters", {})
    mapping.setdefault("filename", "{id}.tif")
    mapping.setdefault("columns", {})
    return mapping


def read_layers(pipeline: str, root: str = ".") -> pd.DataFrame:
    return pd.read_csv(pipeline_dir(pipeline, root) / "layers.csv", keep_default_na=False)


def layer_row(layers: pd.DataFrame, key: str) -> dict:
    rows = layers[layers.filename == f"{key}.tif"]
    if len(rows) != 1:
        raise MappingError(f"Expected one layers.csv row for {key}.tif, found {len(rows)}")
    return rows.iloc[0].to_dict()


# -- value handling -----------------------------------------------------------


def _blank(value) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value)) or value == ""


def _norm(value) -> Optional[str]:
    """Canonical string for comparing layers.csv values with STAC properties:
    ``100``, ``100.0`` and ``"100"`` compare equal; None/NaN/"" are null."""
    if _blank(value):
        return None
    if isinstance(value, bool):
        return str(value).lower()
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if math.isfinite(number) and number == int(number):
        return str(int(number))
    return repr(number)


def _lookup(table: dict, value):
    """Translate ``value`` with a value map whose keys may be strings or numbers."""
    wanted = _norm(value)
    for k, v in (table or {}).items():
        if _norm(k) == wanted or (k is None and wanted is None):
            return v
    return value


def parse_filter(text) -> dict:
    """Parse a ``datapkg_filter`` cell: ``key=value;key=v1|v2``."""
    if _blank(text):
        return {}
    out = {}
    for part in str(text).split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise MappingError(f"Bad datapkg_filter term {part!r}: expected key=value")
        key, value = (s.strip() for s in part.split("=", 1))
        options = [NULL if v == "null" else v for v in value.split("|")]
        out[key] = options if len(options) > 1 else options[0]
    return out


def format_filter(filters: dict) -> str:
    def fmt(v):
        return "null" if v is None else str(v)

    return ";".join(
        f"{k}={'|'.join(map(fmt, v)) if isinstance(v, (list, tuple)) else fmt(v)}"
        for k, v in sorted(filters.items())
    )


def row_filters(row: dict, mapping: dict) -> dict:
    """Property filter selecting the catalogue item for one layers.csv row."""
    filters = dict(mapping.get("filters") or {})
    for column, prop in mapping["properties"].items():
        if column not in row:
            raise MappingError(f"layers.csv has no column {column!r} named in datapkg.yml")
        value = _lookup(mapping["values"].get(column), row[column])
        filters[prop] = None if _blank(value) else value
    filters.update(parse_filter(row.get(FILTER_COLUMN)))
    return filters


def _forms(value) -> list:
    """Accepted string forms of one filter value for the client's string
    comparison: ``100``, ``100.0`` and ``"100"`` all match a property of
    100 or 100.0. None stays None (the client's ``Catalog.ABSENT``)."""
    norm = _norm(value)
    if norm is None:
        return [None]
    forms = {norm, str(value)}
    try:
        number = float(norm)
        if math.isfinite(number):
            forms.add(str(number))
    except ValueError:
        pass
    if norm in ("true", "false"):
        forms.add(norm.capitalize())
    return sorted(forms)


def client_filters(filters: dict) -> dict:
    """``filters`` in the form ``Catalog.items(**filters)`` takes: each value a
    list of accepted strings, None meaning absent or null."""
    out = {}
    for key, wanted in filters.items():
        options = wanted if isinstance(wanted, (list, tuple, set)) else [wanted]
        out[key] = sorted({f for o in options for f in _forms(o)}, key=lambda f: (f is not None, f or ""))
    return out


def matches(properties: dict, filters: dict) -> bool:
    """Local twin of the client's matching (``Catalog.items``) with layers.csv
    number normalisation, for item lists already in hand (refresh_layers).
    None among the accepted values matches an absent or null property."""
    for key, wanted in filters.items():
        options = wanted if isinstance(wanted, (list, tuple, set)) else [wanted]
        actual = _norm(properties.get(key))
        if actual not in {_norm(o) for o in options}:
            return False
    return True


def select(items: list, filters: dict) -> list:
    return [i for i in items if matches(i.properties, filters)]


# -- resolution ------------------------------------------------------------------


def item_spec(row: dict, mapping: dict) -> dict:
    """How a row is resolved, as plain data (used as a Snakemake param, so a
    mapping change re-triggers the fetch job)."""
    dataset = row.get(DATASET_COLUMN)
    dataset = mapping["dataset"] if _blank(dataset) else str(dataset)
    if not _blank(row.get(ITEM_COLUMN)):
        return {"dataset": dataset, "item": str(row[ITEM_COLUMN])}
    return {"dataset": dataset, "filters": row_filters(row, mapping)}


def resolve(catalog, spec: dict):
    """The single catalogue item described by ``spec`` (see :func:`item_spec`)."""
    dataset = spec["dataset"]
    if "item" in spec:
        return catalog.item(dataset, spec["item"])
    # the client filters (and, from a lock, skips reading items its recorded
    # properties rule out); None in a filter matches an absent/null property
    found = catalog.items(dataset, **client_filters(spec["filters"]))
    if len(found) != 1:
        ids = ", ".join(i.id for i in found[:5])
        raise MappingError(
            f"{dataset}: filter {format_filter(spec['filters'])} matched "
            f"{len(found)} items{' (' + ids + ')' if ids else ''}; expected exactly one. "
            f"Set {ITEM_COLUMN} or {FILTER_COLUMN} in layers.csv."
        )
    return found[0]


def locked_sha256(lock: Optional[dict], spec: dict) -> Optional[str]:
    """Checksum the lock file records for a spec's item, if it names an item."""
    if not lock or "item" not in spec:
        return None
    entry = lock.get("datasets", {}).get(spec["dataset"], {}).get("items", {}).get(spec["item"])
    if not entry:
        return None
    data = entry.get("data") or next(iter(entry.values()), {})
    return data.get("sha256")


def read_lock(path) -> Optional[dict]:
    if path and os.path.exists(path):
        with open(path) as f:
            return yaml.safe_load(f)
    return None


def open_catalog(base=None, release=None, lock=None, private_base=None, cache_dir=None):
    """Open the catalogue from a lock file if it exists, else from base/release.

    With a lock file, ``base`` (if given) overrides the base recorded in the
    lock: a mirror of the same catalogue (versions, item ids and checksums
    still come from the lock). ``release`` is then ignored.
    """
    from irv_datapkg.catalog import Catalog

    kwargs = {}
    if cache_dir:
        kwargs["cache_dir"] = cache_dir
    if private_base:
        kwargs["private_base"] = private_base
    if lock and os.path.exists(lock):
        return Catalog.from_lock(lock, base=base or None, **kwargs)
    if not base:
        raise MappingError("datapkg: set 'base' (or commit a lock file) in config.yml")
    return Catalog.open(base, release=release, **kwargs)


def lock_selection(specs) -> tuple:
    """``(datasets, items)`` for ``irv-datapkg lock``: datasets in which some
    layer is resolved by a filter are locked whole (positional DATASETS);
    datasets whose layers all name an item lock only those
    (``--item DATASET/ID``)."""
    by_dataset = {}
    for spec in specs:
        by_dataset.setdefault(spec["dataset"], []).append(spec)
    whole, items = [], []
    for dataset, group in sorted(by_dataset.items()):
        if all("item" in s for s in group):
            items += sorted({f"{dataset}/{s['item']}" for s in group})
        else:
            whole.append(dataset)
    return whole, items


# -- layers.csv generation ------------------------------------------------------


def _reverse(table: dict, value):
    for k, v in (table or {}).items():
        if _norm(v) == _norm(value):
            return "None" if k is None else k
    return value


def item_row(item, mapping: dict) -> dict:
    """Layers-side values of the mapped columns for one catalogue item."""
    row = {}
    for column, prop in mapping["properties"].items():
        value = _reverse(mapping["values"].get(column), item.properties.get(prop))
        row[column] = "" if value is None else value
    return row


class _Default(dict):
    def __missing__(self, key):
        raise MappingError(f"Template field {{{key}}} is not a column, property or 'id'")


def _render(template: str, context: dict) -> str:
    return string.Formatter().vformat(template, (), _Default(context))


def refresh_layers(layers: pd.DataFrame, items: list, mapping: dict, add_missing=False):
    """Fill ``datapkg_item`` in ``layers`` from catalogue items.

    Existing rows keep their ``filename`` (the raster key used downstream) and
    every other column; each is matched to exactly one item. With
    ``add_missing``, items matching no row are appended as new rows, with
    ``filename``/``columns`` rendered from the mapping templates.

    Returns ``(layers, report)`` where report lists unmatched/ambiguous rows
    and unused items.
    """
    layers = layers.copy()
    if ITEM_COLUMN not in layers.columns:
        layers[ITEM_COLUMN] = ""
    report = {"unmatched": [], "ambiguous": [], "new": [], "unused": []}
    used = set()
    for idx, row in layers.iterrows():
        row = row.to_dict()
        if not _blank(row.get(DATASET_COLUMN)) and row[DATASET_COLUMN] != mapping["dataset"]:
            continue  # row from another catalogue dataset: left as is
        if not _blank(row.get(ITEM_COLUMN)):
            ids = [i.id for i in items if i.id == str(row[ITEM_COLUMN])]
        else:
            ids = [i.id for i in select(items, row_filters(row, mapping))]
        if len(ids) == 1:
            layers.at[idx, ITEM_COLUMN] = ids[0]
            used.add(ids[0])
        elif not ids:
            report["unmatched"].append(row["filename"])
        else:
            report["ambiguous"].append((row["filename"], ids))
    new_rows = []
    for item in items:
        if item.id in used:
            continue
        if not add_missing:
            report["unused"].append(item.id)
            continue
        values = item_row(item, mapping)
        context = {**{k: v for k, v in item.properties.items()}, **values, "id": item.id}
        filename = _render(mapping["filename"], context)
        context["stem"] = filename[:-4] if filename.endswith(".tif") else filename
        row = {**values, "filename": filename, ITEM_COLUMN: item.id}
        for column, template in mapping["columns"].items():
            row[column] = _render(str(template), context)
        new_rows.append(row)
        report["new"].append(filename)
    if new_rows:
        layers = pd.concat([layers, pd.DataFrame(new_rows)], ignore_index=True)
        layers = layers.fillna("")
    return layers, report


# -- metadata.json ----------------------------------------------------------------

LOCAL_METADATA_KEYS = ("domain", "name", "group", "keys")


def citation_text(collection: dict) -> Optional[str]:
    citation = collection.get("sci:citation")
    doi = collection.get("sci:doi")
    if citation and doi and doi not in citation:
        citation = f"{citation} https://doi.org/{doi}"
    return citation


def metadata_from_collection(metadata: dict, collection: dict, description="citation") -> dict:
    """Update ``license`` and ``description`` of a pipeline metadata.json from a
    STAC collection; ``domain``, ``name``, ``group`` and ``keys`` stay local.

    The ``raster_tile_sources`` table has no citation column and the frontend
    shows ``description`` as the source credit, so by default ``description``
    takes the collection's ``sci:citation`` (plus DOI); pass
    ``description="description"`` to use the collection description instead.
    """
    out = dict(metadata)
    if collection.get("license"):
        out["license"] = collection["license"]
    text = citation_text(collection) if description == "citation" else None
    text = text or collection.get("description")
    if text:
        out["description"] = text
    for key in LOCAL_METADATA_KEYS:
        if key in metadata:
            out[key] = metadata[key]
    return out


def dump_metadata(metadata: dict) -> str:
    return json.dumps(metadata, indent=2, ensure_ascii=False) + "\n"
