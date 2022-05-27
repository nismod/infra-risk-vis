FROM python:3.10
LABEL maintainer="frederick.thomas@ouce.ox.ac.uk"

WORKDIR /code

COPY Pipfile* .

COPY . /code
RUN pip install pipenv
RUN pipenv install

CMD ["pipenv", "run", "uvicorn", "backend.app.main:app", "--host", "0.0.0.0", "--port", "8888"]

# If running behind a proxy like Nginx or Traefik add --proxy-headers
# CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80", "--proxy-headers"]
