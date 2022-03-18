FROM python:3.10

WORKDIR /code

RUN pip install --no-cache-dir --upgrade pipenv

COPY ./Pipfile /code/Pipfile
COPY ./Pipfile.lock /code/Pipfile.lock

RUN pipenv install --deploy --ignore-pipfile

COPY ./backend /code/backend

CMD ["pipenv", "run", "uvicorn", "backend.app.main:app", "--host", "0.0.0.0", "--port", "80"]

# If running behind a proxy like Nginx or Traefik add --proxy-headers
# CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80", "--proxy-headers"]