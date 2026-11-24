# Proposed GitHub Actions workflows

Workflows that could not be committed under `.github/workflows/` from the
session that wrote them. This repository has no workflows yet. A maintainer
should copy them in:

```bash
mkdir -p .github/workflows
cp ci/etl-datapkg.yml .github/workflows/etl-datapkg.yml
git rm -r ci/
```

| File | What it runs |
| --- | --- |
| `etl-datapkg.yml` | `etl/tests` (irv-datapkg integration of the raster ETL) on Python 3.13 with pip wheels (`snakemake pytest pandas xarray` + irv-datapkg from git at `IRV_DATAPKG_REF`, default `main`; a workflow_dispatch input overrides it). Runs on changes under `etl/` and weekly. No database or terracotta service is needed. |

Until irv-datapkg 0.3.0 is released, run it with `irv_datapkg_ref` set to the
irv-datapkg integration branch (or change the default in `env:`).
