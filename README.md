# Infrastructure Risk Visualisation Tool

This project provides interactive data visualisations of risk analysis results.

![About](images/screenshot-about.png)

The tool presents the infrastructure systems and hazards considered in the
analysis, then presents results as modelled for the whole system at a fine
scale.

See an overview of infrastructure networks, hazards, risks:

![Networks](images/screenshot-overview.png)

Other core functionality:

- Click road or rail links to see attributes, assessed risks of direct damages
  to assets and indirect economic losses resulting from infrastructure service
  disruption;
- Explore a cost-benefit analysis of adaptation measures;
- Assess the wider sustainability of adaptation and mitigation options.

This README covers requirements and steps through how to prepare data for
visualisation and how to run the tool.

1. Data preparation
2. Build and run
3. Deployment

## Data preparation

The visualisation tool runs using prepared versions of analysis data and results

- Rasters stored as Cloud-Optimised GeoTIFFs, with metadata ingested into
  a terracotta SQLite database
- Vector data stored in a PostgreSQL database, and preprocessed into Mapbox
  Vector Tiles

See `./etl` directory for details.

## Build and run

Running the application requires several (local) server processes: the
vector and raster tileservers, the app backend, and the app frontend.

### Node and npm

The build and run steps use [node.js](https://nodejs.org/) - this provides the
`npm` command.

Install required packages. Run from the project root:

    npm install

### Install the raster tileserver

Install the raster tileserver - [`terracotta`](https://terracotta-python.readthedocs.io/en/latest/)

For example, installing using conda:

    conda create --name infrariskvis python=3.8 numpy rasterio shapely crick
    conda activate infrariskvis
    pip install terracotta[recommended]

### Run the vector tileserver

Install the vector tileserver - [`tileserver-gl-light`](https://tileserver.readthedocs.io/en/latest/)

    npm install tileserver-gl-light

Run the vector tileserver (from the root of the project):

    npm run vector

### Run the raster tileserver

Prepare the raster tileserver database, providing the path to the directory that
contains the hazard GeoTIFFs:

    npm run raster-init ./path/to/raster/data/

Run the raster tileserver:

    npm run raster

### Run the backend API server and database

Two options here, with or without docker.

Without docker, follow the notes in `./backend/README.md` to setup a development
environment for python.

Set up a postgres database and add connection details in `./backend/.env`.

Run the api server:

    cd ./backend
    pipenv run uvicorn backend.app.main:app --host localhost --port 8888

Alternatively, run `docker compose up` to run the API server in one container and
postgres in another.

### Run the frontend app in development mode

Start the app server:

    npm start

This should automatically open a browser tab. If not, open:

    firefox http://localhost:3000/

## Deployment

The site can run on a single Linux machine or virtual machine, with a suggested
configuration that deploys the server processes behind an Nginx reverse proxy
in production modes.

See `./deploy` directory for details.

## Acknowledgements

This tool has been developed through several projects.

- [v0.1](https://github.com/oi-analytics/oi-risk-vis/releases/tag/v0.1-argentina)
  was developed by Oxford Infrastructure Analytics for the Government of
  Argentina with funding support from the World Bank Group and Global Facility
  for Disaster Reduction and Recovery (GFDRR).
- [v0.2](https://github.com/oi-analytics/oi-risk-vis/releases/tag/v0.2.0-seasia)
  was developed by Oxford Infrastructure Analytics for the Disaster Risk
  Financing and Insurance Program (DRFIP) of the World Bank with support from
  the Japan&mdash;World Bank Program for Mainstreaming DRM in Developing
  Countries, which is financed by the Government of Japan and managed by the
  Global Facility for Disaster Reduction and Recovery (GFDRR) through the Tokyo
  Disaster Risk Management Hub.
- [v0.3](https://github.com/oi-analytics/oi-risk-vis/releases/tag/v0.3.0-jamaica)
  was developed by the Oxford Programme for Sustainable Infrastructure Systems
  (OPSIS) in the Environmental Change Institute, University of Oxford, for the
  Government of Jamaica (GoJ) as part of a project funded by UK Aid (FCDO). The
  initiative forms part of the Coalition for Climate Resilient Investment’s
  (CCRI) collaboration with the GoJ, which also includes analysis of
  nature-based approaches to build resilience in Jamaica to be procured and
  funded by the Green Climate Fund (GCF).
- [release/caribbean](https://github.com/nismod/infra-risk-vis/tree/release/caribbean)
  was developed as part of the Jamaica project.
- [release/east-africa](https://github.com/nismod/infra-risk-vis/tree/release/east-africa)
  was developed by researchers in the University of Southampton's Transportation
  Research Group and the Oxford Programme for Sustainable Infrastructure
  Systems, University of Oxford, supported by engagement with infrastructure and
  climate specialists and related government bodies, and funded by UKAID through
  the UK Foreign, Commonwealth & Development Office under the High Volume
  Transport Applied Research Programme, managed by DT Global.
