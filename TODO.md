- TLS support
    - certbot
    - Self-signed certs for development? Necessary?
- ENV file
    - Propagate DB conn. parameters into BE containers
    - Can we pull out any more config?
- Development stack
    - FE code mount
    - FE hot reload?
    - BE Pipfile development packages?
    - BE code mount
    - Move DB service into compose-local.yml
- Host provisioning requirements
    - docker, docker compose
- Slow to die containers
    - database
    - raster-tileserver
- Stop using virtual environment tooling
    - pipenv in the backend
    - conda in raster-tileserver, see following from deployment script:
        sudo python3.10 -m pip install --upgrade cython
        sudo python3.10 -m pip install --upgrade numpy
        sudo python3.10 -m pip install --upgrade --no-binary rasterio rasterio==1.3a4
        sudo python3.10 -m pip install --upgrade gunicorn terracotta[recommended]