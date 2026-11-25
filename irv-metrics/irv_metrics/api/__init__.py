"""The /v1 metrics API (contract C1): `Service` (reads, no HTTP) and `create_router` (FastAPI)."""
from .routes import create_router
from .service import ApiError, Service

__all__ = ["ApiError", "Service", "create_router"]
