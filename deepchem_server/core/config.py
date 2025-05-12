from typing import Union
from deepchem_server.core.datastore import DiskDataStore

_DATASTORE = None


def set_datastore(datastore):
    global _DATASTORE
    _DATASTORE = datastore


def get_datastore() -> Union[DiskDataStore, None]:
    global _DATASTORE
    return _DATASTORE


def refresh():
    set_datastore(None)
