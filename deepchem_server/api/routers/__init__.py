"""API Routers package.

Contains all versioned API routers.
"""
from deepchem_server.api.routers.v1 import router as v1_router

__all__ = ["v1_router"]
