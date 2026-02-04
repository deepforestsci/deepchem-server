"""
Datastore Client Package - Client library for accessing the datastore service.

Provides:
- DatastoreClient: Low-level HTTP client for REST API
- DeepchemDatastore: High-level interface used by primitives

Usage:
    from deepchem_server.services.datastore.client import (
        DatastoreClient,
        DeepchemDatastore,
    )
    
    # Create client
    client = DatastoreClient(url="http://localhost:8081")
    
    # Create datastore instance for a profile/project
    datastore = DeepchemDatastore(client, profile="user", project="my_project")
    
    # Use like before
    datastore.upload_data("file.csv", "/path/to/file.csv", card)
    data = datastore.get("deepchem://user/my_project/file.csv")

Future extensibility:
    # AWS S3 backend (future)
    client = S3DatastoreClient(bucket="my-bucket", prefix="deepchem/")
    datastore = DeepchemDatastore(client, profile="user", project="my_project")
"""
from deepchem_server.services.datastore.client.datastore import DeepchemDatastore
from deepchem_server.services.datastore.client.http import DatastoreClient


__all__ = ["DatastoreClient", "DeepchemDatastore"]
