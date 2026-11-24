FROM python:3.13
LABEL maintainer="frederick.thomas@ouce.ox.ac.uk"

WORKDIR /code

# pg_isready, for services that wait for the database (metrics-load)
RUN apt-get update \
    && apt-get install -y --no-install-recommends postgresql-client \
    && rm -rf /var/lib/apt/lists/*

# irv-metrics (the national metrics API, ../irv-metrics) comes from the named
# build context `irv-metrics` (compose: additional_contexts). It is installed
# first, editable, so the backend's plain-name dependency `irv-metrics` is
# already satisfied and pip doesn't go to PyPI (yet). A change to irv-metrics
# rebuilds from here. The dev compose file mounts ./irv-metrics over
# /code/irv-metrics.
COPY --from=irv-metrics . /code/irv-metrics
RUN python -m pip install -e /code/irv-metrics

# Install project with dependencies - these can be cached as docker layers if unchanged
COPY pyproject.toml .
RUN python -m pip install .

# Copy and install latest package code
COPY ./backend /code/backend
RUN python -m pip install --upgrade --compile .

# Database migrations: run with `alembic upgrade head` from /code
COPY alembic.ini /code/alembic.ini
COPY ./migrations /code/migrations

CMD ["uvicorn", "backend.app.main:app", "--host", "0.0.0.0", "--port", "8888"]

# If running behind a proxy like Nginx or Traefik add --proxy-headers
# CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80", "--proxy-headers"]
