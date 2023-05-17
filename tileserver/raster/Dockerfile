FROM python:3.10

WORKDIR /code
RUN python -m pip install --upgrade pip setuptools wheel
RUN pip install cython "numpy>=1.15,!=1.17.0,<1.24"
RUN pip install gunicorn terracotta[recommended]~=0.8

CMD ["gunicorn", "--workers", "3", "--bind", "0.0.0.0:5000", "-m", "007", "--access-logfile", "'-'", "terracotta.server.app:app"]
