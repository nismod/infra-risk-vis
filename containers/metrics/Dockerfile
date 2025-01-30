FROM python:3.13
LABEL maintainer="alexander.shiarella@humanities.ox.ac.uk"

WORKDIR /code

# Install project with dependencies - these can be cached as docker layers if unchanged
COPY pyproject.toml .
RUN python -m pip install .

# Copy and install latest package code
COPY ./api /code/api

RUN python -m pip install --upgrade --compile .

CMD ["uvicorn", "api.app.main:app", "--host", "0.0.0.0", "--port", "8888"]

# If running behind a proxy like Nginx or Traefik add --proxy-headers
# CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80", "--proxy-headers"]
