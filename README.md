# oi-risk-vis
Risk analysis visualisation tool

https://www.mapbox.com/mapbox-gl-js/style-spec/
Maputnik for developing syles

# Prepare data
Download the ``boundaries`` and ``network`` data

Unzip in ``/incoming_data`` folder

    unzip ~/Downloads/boundaries.zip -d incoming_data/
    unzip ~/Downloads/network.zip -d incoming_data/

Convert the incoming shapefiles into a *.mbtiles file

    ./convert.sh

Download fonts

    wget -P ~/Downloads/ "https://github.com/klokantech/tileserver-gl-styles/archive/master.zip"

# Run the tile server

cd data

docker run -it -v $(pwd):/data -p 8080:80 klokantech/tileserver-gl --config config.json