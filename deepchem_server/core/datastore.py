"""
DEPRECATED: This module is deprecated and will be removed.

The DiskDataStore class has been removed. All storage now goes through
the centralized datastore service.

Migration:
    # Old (removed)
    from deepchem_server.core.datastore import DiskDataStore
    datastore = DiskDataStore(profile, project, basedir)
    
    # New
    from deepchem_server.services.datastore.client import DatastoreClient, DeepchemDatastore
    client = DatastoreClient(url=os.getenv("DATASTORE_URL"))
    datastore = DeepchemDatastore(client, profile, project)

For tests, either:
1. Mock the datastore, or
2. Run a local datastore service
"""
import warnings

warnings.warn(
    "deepchem_server.core.datastore is deprecated and will be removed. "
    "Use deepchem_server.services.datastore.client instead. "
    "See module docstring for migration guide.",
    DeprecationWarning,
    stacklevel=2,
)

# Re-export for minimal backward compatibility during transition
from deepchem_server.services.datastore.client import DeepchemDatastore

# Alias for code that might reference the old name
DiskDataStore = DeepchemDatastore
DataStore = DeepchemDatastore

__all__ = ["DiskDataStore", "DataStore", "DeepchemDatastore"]
