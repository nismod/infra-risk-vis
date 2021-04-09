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

clean:
	rm -f ./data/*.mbtiles
