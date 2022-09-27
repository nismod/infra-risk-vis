FROM python:3.10

RUN apt-get update && \
  apt-get install -y build-essential gdal-bin libgdal-dev && \
  pip install cython && \
  git clone --depth 1 --branch v0.7.5 https://github.com/DHI-GRAS/terracotta && \
  cd terracotta && \
  pip install -e . && \
  pip install gunicorn crick

COPY ./docker-entrypoint.sh /
ENTRYPOINT [ "/docker-entrypoint.sh" ]