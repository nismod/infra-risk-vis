FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
ENV DEBIAN_FRONTEND noninteractive

# Step 1: Retrieve conda environments

# Step 2: Generate conda environments

# RUN apt-get update && apt-get install -y snakemake && \
RUN mamba clean --all -y && \
    mamba create -y -c bioconda -c conda-forge -n snakemake snakemake-minimal
# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "snakemake", "/bin/bash", "-c"]
RUN mamba install -y geopandas && \
    mamba install -y rasterio && \
    mamba install -y gdal && \
    mamba install -y terracotta && \
    mamba install -y tippecanoe && \
    mamba install -y jq && \
    mamba install -y pyarrow && \
    mamba install -y psycopg2-binary sqlalchemy==1.4.41 geoalchemy2==0.12.5 ujson==5.5.0

# Test the environment is ready:
RUN echo "Testing environment..."
RUN python -c "import geopandas" > /dev/null

# Ensure the ETL folder is setup and the working dir
RUN mkdir /opt/etl
WORKDIR /opt/etl

# Setup default image env
ENV CONDA_DEFAULT_ENV snakemake

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "snakemake", "snakemake"]

