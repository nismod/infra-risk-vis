FROM python:3.12
LABEL maintainer="frederick.thomas@ouce.ox.ac.uk"

WORKDIR /code

# Install project with dependencies - these can be cached as docker layers if unchanged
COPY pyproject.toml .
RUN python -m pip install .

# Copy and install latest package code
COPY ./backend /code/backend
RUN python -m pip install --upgrade --compile .

CMD ["uvicorn", "backend.app.main:app", "--host", "0.0.0.0", "--port", "8888"]

# If running behind a proxy like Nginx or Traefik add --proxy-headers
# CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80", "--proxy-headers"]
