FROM condaforge/mambaforge:latest
ENV DEBIAN_FRONTEND noninteractive

# Generate conda environments

RUN mamba clean --all -y && \
    mamba create -y -c bioconda -c conda-forge -n snakemake snakemake-minimal
# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "snakemake", "/bin/bash", "-c"]
RUN mamba install -y gdal geopandas rasterio rioxarray xarray netcdf4 terracotta tippecanoe jq pyarrow \
    psycopg2-binary sqlalchemy~=1.4 geoalchemy2 ujson

# Test the environment is ready:
RUN echo "Testing environment..."
RUN python -c "import geopandas" > /dev/null

# Ensure the ETL folder is setup and the working dir
RUN mkdir /opt/etl
WORKDIR /opt/etl

# Setup default image env
ENV CONDA_DEFAULT_ENV snakemake

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "snakemake", "snakemake"]
