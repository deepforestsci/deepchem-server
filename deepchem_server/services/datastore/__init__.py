"""
DeepChem Datastore Package

Centralized file storage solution with HTTP-based access.

Sub-packages:
    server - DatastoreService, FastAPI REST API (run as standalone service)
    client - DatastoreClient, DeepchemDatastore (use in workers/gateway)

Usage:
    # Server (run as service)
    python -m deepchem_server.services.datastore
    
    # Client (in code)
    from deepchem_server.services.datastore.client import DatastoreClient, DeepchemDatastore
    
    client = DatastoreClient(url="http://localhost:8081")
    datastore = DeepchemDatastore(client, profile="user", project="project")
"""

# Re-export for convenience
from deepchem_server.services.datastore.client import DatastoreClient, DeepchemDatastore
from deepchem_server.services.datastore.server import DatastoreService


__all__ = [
    "DatastoreService",
    "DatastoreClient",
    "DeepchemDatastore",
]
