# DeepChem Services Package
"""
Services for DeepChem Server.

Sub-packages:
    datastore - Centralized file storage service (server and client)
    jobs - Job queue and worker management

Usage:
    # Run datastore service
    python -m deepchem_server.services.datastore
    
    # Import client
    from deepchem_server.services.datastore import DatastoreClient, DeepchemDatastore
"""
from deepchem_server.services.datastore import (
    DatastoreClient,
    DatastoreService,
    DeepchemDatastore,
)


__all__ = ["DatastoreService", "DatastoreClient", "DeepchemDatastore"]
