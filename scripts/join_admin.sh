set -e
set -x

# copy summary across
cp ../seasia/results/expected_summary.csv ./intermediate_data/

# join to existing admin mbtiles
tile-join \
  -pk \
  -o data/admin1_data.mbtiles \
  -c intermediate_data/expected_summary.csv \
  data/admin1.mbtiles

# decode from mbtile
tippecanoe-decode data/admin1_data.mbtiles 3 6 4 > data/admin1_364.json
tippecanoe-decode data/admin1_data.mbtiles 3 6 3 > data/admin1_363.json
tippecanoe-decode data/admin1_data.mbtiles 3 7 4 > data/admin1_374.json
tippecanoe-decode data/admin1_data.mbtiles 3 7 3 > data/admin1_373.json

# concatenate
cat data/admin1_*.json | grep '"Feature"' > data/admin1_data.json

# use QGIS - merge on GID_1

# rebuild mbtiles with single layer
tippecanoe \
  -zg \
  --output=./data/admin1_merged.mbtiles \
  --layer=admin1 \
  ./data/admin1_merged.json
