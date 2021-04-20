.PHONY: all clean

all: ./data/rail.mbtiles ./data/roads_main.mbtiles ./data/roads_other.mbtiles ./data/electricity.mbtiles

./data/rail.mbtiles:
	tippecanoe \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--output=./data/rail.mbtiles \
		--layer=rail \
		./intermediate_data/Rail/*.json

./data/roads_main.mbtiles:
	tippecanoe \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-smallest-as-needed \
		--output=./data/roads_main.mbtiles \
		./intermediate_data/Roads/highway/trunk.json \
		./intermediate_data/Roads/highway/motorway.json \
		./intermediate_data/Roads/highway/primary.json \
		./intermediate_data/Roads/highway/secondary.json

./data/roads_other.mbtiles:
	tippecanoe \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-densest-as-needed \
		--output=./data/roads_other.mbtiles \
		--layer=other \
		--read-parallel \
		./intermediate_data/Roads/highway/tertiary.json \
		./intermediate_data/Roads/highway/other_*.json

./data/electricity.mbtiles:
	tippecanoe \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-smallest-as-needed \
		--output=./data/electricity.mbtiles \
		--layer=electricity \
		./intermediate_data/Electricity/*.geobuf

./data/coastal.mbtiles:
	tippecanoe \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-rate=1 \
		--cluster-distance=1 \
		--accumulate-attribute=depth_m:max \
		--output=./data/coastal.mbtiles \
		--layer=coastal \
		./intermediate_data/Flooding/coastal_100yr.csv

./data/fluvial.mbtiles:
	tippecanoe \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-rate=1 \
		--cluster-distance=1 \
		--accumulate-attribute=depth_m:max \
		--output=./data/fluvial.mbtiles \
		--layer=fluvial \
		./intermediate_data/Flooding/fluvial_100yr.csv

./data/cyclone.mbtiles:
	tippecanoe \
		-zg \
		--output=./data/cyclone.mbtiles \
		--layer=cyclone \
		./intermediate_data/Cyclone/cyclone_100yr.json

clean:
	rm -f ./data/*.mbtiles

./intermediate_data/Cyclone/cyclone_100yr.json:
	gdal_polygonize.py Cyclone_100yr.tif -f "GeoJSONSeq" cyclone_100yr.json cyclone gust_speed_ms-1
