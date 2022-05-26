from uvicorn.workers import UvicornWorker

class ProductionUvicornWorker(UvicornWorker):
    CONFIG_KWARGS = {
        "forwarded_allow_ips": '*',
        "proxy_headers": True,
        "uds": "backend.sock"
    }
