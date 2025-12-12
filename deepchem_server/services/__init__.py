# DeepchemDatastore Standalone Service
"""
Standalone HTTP/REST service for file storage operations.

Usage:
    python -m deepchem_server.services.datastore
"""
from deepchem_server.services.datastore.service import DatastoreService
from deepchem_server.services.datastore.client import DatastoreClient

__all__ = ["DatastoreService", "DatastoreClient"]
