#
# Convert flood tiffs to CSV for mbtiles generation
#
pushd intermediate_data/Flooding
  # convert to points space-separated text file (x, y, value columns, no header)
  gdal_translate -of XYZ Coastal_100yr.tif coastal_100yr.txt
  gdal_translate -of XYZ Fluvial_100yr.tif fluvial_100yr.txt

  # drop zero values
  grep -v " 0$" fluvial_100yr.txt > fluvial_100yr_nonzero.txt
  grep -v " 0$" coastal_100yr.txt > coastal_100yr_nonzero.txt

  # drop large files
  rm *yr.txt

  # rename to CSV
  rename 's/_nonzero.txt//' *.csv

  # swap spaces for commas
  sed -i 's/ /,/g' *.csv

  # add header lines
  echo 'x,y,depth_m' | cat - fluvial_100yr.csv > _fluvial_100yr.csv \
    && mv _fluvial_100yr.csv fluvial_100yr.csv

  echo 'x,y,depth_m' | cat - coastal_100yr.csv > _coastal_100yr.csv \
    && mv _coastal_100yr.csv coastal_100yr.csv

popd
