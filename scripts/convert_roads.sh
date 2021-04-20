#
# Convert JSON Roads dump to files per road type, for useful layers as mbtiles
#
set -e
set -x

pushd intermediate_data/Roads

  # take country JSON files and filter to file-per-road-type
  # grep for "primary" etc. will include "primary_link", which we *do* want
  grep -h trunk *.json > highway/trunk.json
  grep -h motorway *.json > highway/motorway.json
  grep -h primary *.json > highway/primary.json
  grep -h secondary *.json > highway/secondary.json
  grep -h tertiary *.json > highway/tertiary.json

  # everything that's non-main, which is a lot!
  for f in *.json
  do
    grep highway $f | egrep -v "trunk|motorway|primary|secondary|tertiary" > highway/other_$f
  done

popd
