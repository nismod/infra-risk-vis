# OIA Risk Visualisation Tool

This project provides interactive data visualisations of OIA risk analysis results.

This README covers requirements and steps through how to prepare data for visualisation and how
to run the tool.

1. Data preparation requirements
2. Prepare data
3. Build and run requirements
4. Run


## Data preparation requirements

### ogr2ogr

[ogr2ogr](https://www.gdal.org/ogr2ogr.html) is used for spatial data processing. On Ubuntu,
run:

    sudo apt-get install gdal-bin

### Tippecanoe

The data preparation steps use [Mapbox tippecanoe](https://github.com/mapbox/tippecanoe) to
build vector tiles from large feature sets.

The easiest way to install tippecanoe on OSX is with Homebrew:

    brew install tippecanoe

On Ubuntu it will usually be easiest to build from the source repository:

    sudo apt-get install build-essential g++ libsqlite3-dev zlib1g-dev
    git clone https://github.com/mapbox/tippecanoe
    cd tippecanoe
    make -j
    make


## Prepare data

For Argentina (for example) download `boundaries`, `network` and `flood_data` from the OIA
shared folder `302 Argentina/D Work Processes/Argentina/data/`.

Unzip within `/incoming_data` folder:

    unzip ~/Downloads/boundaries.zip -d incoming_data/
    unzip ~/Downloads/network.zip -d incoming_data/

Convert the incoming shapefiles into a *.mbtiles file:

    make


## Build and run requirements


### Node and npm

The build and run steps use [node.js](https://nodejs.org/) - this provides the `npm` command.

Install required packages. Run from the project root:

    npm install

### Tileserver GL

Install tlieserver-gl globally:

    npm install -g tileserver-gl

### Docker (optional)

Alternately, Docker CE may be used to run the tileserver in a container. Follow the [docker
installation guide](https://docs.docker.com/install/) to get set up.

## Run

Running the application currently requires two (local) server processes: the tileserver and the
app itself.

### Run the tileserver

Run the tileserver directly (from the root of the project):

    tileserver-gl

Or start the docker container:

    docker run -it -v $(pwd):/data -p 8080:80 klokantech/tileserver-gl --config config.json

Open a browser to view the tileserver:

    firefox http://localhost:8080/

### Run the app

Start the app server:

    npm start

This should automatically open a browser tab. If not, open:

    firefox http://localhost:3000/
