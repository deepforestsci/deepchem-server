from __future__ import annotations

import os
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from deepchem_server.services.datastore.client import DeepchemDatastore


_DATASTORE: Optional["DeepchemDatastore"] = None

# Environment variable for remote datastore URL
DATASTORE_URL = os.getenv("DATASTORE_URL")
DATASTORE_API_KEY = os.getenv("DATASTORE_API_KEY")


def set_datastore(datastore: Optional["DeepchemDatastore"]) -> None:
    """Set the global datastore instance.

    Parameters
    ----------
    datastore : DeepchemDatastore or None
        The datastore instance to set as the global datastore, or None to reset.

    Returns
    -------
    None
    """
    global _DATASTORE
    _DATASTORE = datastore


def get_datastore() -> Optional["DeepchemDatastore"]:
    """Get the current global datastore instance.

    Returns
    -------
    DeepchemDatastore or None
        The current datastore instance, or None if no datastore has been set.
    """
    return _DATASTORE


def get_datastore_client():
    """Get a DatastoreClient if remote mode is configured.
    
    Returns None if DATASTORE_URL is not set.
    
    Returns
    -------
    DatastoreClient or None
        Client for remote datastore service
    """
    if not DATASTORE_URL:
        return None

    from deepchem_server.services.datastore.client import DatastoreClient
    return DatastoreClient(url=DATASTORE_URL, api_key=DATASTORE_API_KEY)


def is_remote_datastore() -> bool:
    """Check if remote datastore mode is enabled."""
    return DATASTORE_URL is not None


def refresh() -> None:
    """Reset the global datastore to None.

    Returns
    -------
    None
    """
    set_datastore(None)
