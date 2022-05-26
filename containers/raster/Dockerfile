FROM continuumio/miniconda3
LABEL maintainer="frederick.thomas@ouce.ox.ac.uk"

RUN mkdir -p /opt/terracotta/db
WORKDIR /opt/terracotta

# install dependencies
#RUN apt-get update && apt-get install -y build-essential gdal-bin libgdal-dev

# it was easiest to get rasterio working (linked with gdal) when using conda
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority strict
COPY environment.yml .
RUN conda env create -f environment.yml

# make RUN commands use the new environment:
#SHELL ["conda", "run", "-n", "raster-tileserver", "/bin/bash", "-c"]

# copy across config and scripts to working directory
COPY config.toml entrypoint.sh .

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "raster-tileserver", "/bin/bash", "entrypoint.sh"]
