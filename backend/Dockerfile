FROM python:3.10

WORKDIR /code

RUN pip install --no-cache-dir --upgrade pipenv

COPY ./Pipfile /code/Pipfile
COPY ./Pipfile.lock /code/Pipfile.lock
COPY ./setup.py /code/setup.py
COPY ./backend /code/backend

RUN pipenv install --deploy --ignore-pipfile

CMD ["pipenv", "run", "uvicorn", "backend.app.main:app", "--host", "0.0.0.0", "--port", "80"]
