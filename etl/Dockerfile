FROM condaforge/mambaforge:latest
ENV DEBIAN_FRONTEND noninteractive

# Generate conda environment
COPY ./environment.yml environment.yml
RUN mamba clean --all -y && \
    mamba create -y -f environment.yml
# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "irv-etl", "/bin/bash", "-c"]

# Ensure the ETL folder is setup and the working dir
RUN mkdir /opt/etl
WORKDIR /opt/etl

# Setup default image env
ENV CONDA_DEFAULT_ENV irv-etl

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "irv-etl", "snakemake"]
