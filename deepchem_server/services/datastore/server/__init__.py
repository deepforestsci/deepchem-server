"""
Datastore Server Package - HTTP/REST service for file storage.

This package provides the centralized file storage service that all
workers and services use via HTTP.

Usage:
    python -m deepchem_server.services.datastore.server
    
Or import for embedding:
    from deepchem_server.services.datastore.server import create_app, DatastoreService
"""
from deepchem_server.services.datastore.server.api import create_app, router
from deepchem_server.services.datastore.server.service import DatastoreService


__all__ = ["DatastoreService", "create_app", "router"]
