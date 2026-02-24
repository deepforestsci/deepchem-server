from __future__ import annotations

from typing import Optional

from deepchem_server.core.datastore import DiskDataStore


_DATASTORE: Optional[DiskDataStore] = None


def set_datastore(datastore: Optional[DiskDataStore]) -> None:
    """Set the global datastore instance.

    Parameters
    ----------
    datastore : DiskDataStore or None
        The datastore instance to set as the global datastore, or None to reset.

    Returns
    -------
    None
    """
    global _DATASTORE
    _DATASTORE = datastore


def get_datastore() -> Optional[DiskDataStore]:
    """Get the current global datastore instance.

    Returns
    -------
    DiskDataStore or None
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
