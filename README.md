# Infrastructure Risk Visualisation Tool

This project provides interactive data visualisations of risk analysis results.

![About](images/screenshot-about.png)

The tool presents the infrastructure systems and hazards considered in the
analysis, then presents results as modelled for the whole system at a fine
scale.

See an overview of infrastructure networks:

![Networks](images/screenshot-overview.png)

Other functionality planned (and incorporated in some way in previous versions):

- Summarise risk analysis at an administrative regional scale.
- Zoom in to see networks in detail.
- See an overview of hazard data.
- Inspect details of hazard layers.
- Query attributes of elements of the system.
- Range of potential economic impacts of failure, consisting of direct damages
  to infrastructure assets and indirect economic losses resulting from
  infrastructure service disruption (loss of power, loss of access).
- Explore a cost-benefit analysis (under uncertainty, with options to explore
  some parameters) of adaptation measures.

This README covers requirements and steps through how to prepare data for
visualisation and how to run the tool.

1. Data preparation
3. Build and run
4. Deployment

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

### Terracotta

Install the raster tileserver - Terracotta

For example, installing using conda:

    conda create --name infrariskvis python=3.8 numpy rasterio shapely crick
    conda activate infrariskvis
    pip install terracotta[recommended]

### Run the vector tileserver

Run the tileserver directly (from the root of the project):

    npm run vector

### Run the raster tileserver

Prepare the raster tileserver database:

    npm run raster-init

Run the raster tileserver:

    npm run raster

### Run the backend API server and database

Two options here.

Without docker, follow the notes in `./backend/README.md` to setup a development
environment for python.

Set up a postgres database and add connection details in `./backend/.env`.

Run the api server:

    cd ./backend
    pipenv run uvicorn backend.app.main:app --host localhost --port 8888

Alternatively, run `docker-compose` to run the API server in one container and
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
- current development is by the Oxford Programme for Sustainable Infrastructure
  Systems in the Environmental Change Institute, University of Oxford, for the
  Government of Jamaica (GoJ) as part of a project funded by UK Aid (FCDO). The
  initiative forms part of the Coalition for Climate Resilient Investment’s
  (CCRI) collaboration with the GoJ, which also includes analysis of
  nature-based approaches to build resilience in Jamaica to be procured and
  funded by the Green Climate Fund (GCF).
