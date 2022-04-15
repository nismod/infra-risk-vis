# API backend

Server app, written in Python, includes database definition and etl.

`Dockerfile` shows example setup. See
[docs](https://pipenv.pypa.io/en/latest/basics/#pipenv-and-docker-containers)
for how to build venv and run in a more production-oriented way.

Outline of dependencies:
- system `python` and `pip` (e.g. `apt install python3 python3-pip` on Ubuntu)
- `pipenv` to manage dependencies (e.g. `pip install pipenv`)
- optionally [`pyenv`](https://github.com/pyenv/pyenv) to provide alternative
  version to system python
- `Pipfile` defines python package dependencies
  - includes editable install of `backend` for this application as a package, as
    minimally configured in `setup.py`
- system library dependencies of python packages include
  `libgdal-dev libgeos-dev libpq-dev libproj-dev`

Using `pipenv`:
- run `pipenv install --dev` to install into dev environment
- run `pipenv run ...` to run individual commands within the environment:
  - `pipenv run python` for a REPL
  - `pipenv run snakemake -s ../etl/Snakefile --cores 1` for snakemake
  - `pipenv run uvicorn backend.app.main:app --host locahost --port 8888`
    for the API server
- run `pipenv shell` to drop into a shell within the environment

Environment variables:
- can use `.env` to define environment variables, `pipenv` will load them
  automatically. Use
  [`PG*`](https://www.postgresql.org/docs/current/libpq-envars.html) to define
  database connection details. See `.env.example` for an example.
