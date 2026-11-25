Fixture for `tests/test_csv_loader.py`: a dataset config and a data directory
for `irv_metrics.loaders.csv`.

- `config/`: a vocabulary with one admin set (`test-admin`, two GeoJSON level
  files) and two datasets (with ids from the GRI config; the loader seeds the
  cube from this config).
- `data/ZZA`: an admin-0 file whose `adm_id` is the country name, an admin-1
  file with one empty value, and a concentration curve (wide, 3 points).
- `data/ZZB`: a blank admin-0 stub next to an admin-1 file, which the admin-0
  completeness rule rejects.
- `data/boundaries`: ZZA (two admin-1 areas), ZZB, and in the admin-1 file a
  coarser `shapeType: ADM0` polygon that must not be loaded.
