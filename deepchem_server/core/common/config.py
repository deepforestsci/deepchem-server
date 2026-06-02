from __future__ import annotations

from typing import Optional

from deepchem_server.core.datastore import DataStore


_DATASTORE: Optional[DataStore] = None


def set_datastore(datastore: Optional[DataStore]) -> None:
    """Set the global datastore instance.

    Parameters
    ----------
    datastore : DataStore or None
        The datastore instance to set as the global datastore, or None to reset.

    Returns
    -------
    None
    """
    global _DATASTORE
    _DATASTORE = datastore


def get_datastore() -> Optional[DataStore]:
    """Get the current global datastore instance.

    Returns
    -------
    DataStore or None
        The current datastore instance, or None if no datastore has been set.
    """
    return _DATASTORE


def refresh() -> None:
    """Reset the global datastore to None.

    Returns
    -------
    None
    """
    set_datastore(None)
