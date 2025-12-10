import multiprocessing as mp

import pandas as pd
import pytest

from deepchem_server.core.datastore import DiskDataStore


@pytest.fixture(scope="session", autouse=True)
def set_multiprocessing_start_method():
    """Set multiprocessing start method to 'fork' for the test session.

    On macOS, the default is 'spawn' which creates fresh Python processes
    that don't inherit global state. This causes 'Datastore not set' errors
    in multicore featurization tests. Using 'fork' ensures child processes
    inherit the global datastore configuration.
    """
    try:
        mp.set_start_method('fork', force=True)
    except RuntimeError:
        # Already set, ignore
        pass


@pytest.fixture
def disk_datastore(tmp_path):
    fds = DiskDataStore(profile_name='test', project_name='user', basedir=str(tmp_path))
    return fds


@pytest.fixture
def alternate_disk_datastore(tmp_path):
    fds = DiskDataStore(profile_name='alternate-test', project_name='alternate-user', basedir=str(tmp_path))
    return fds


@pytest.fixture
def tmp_csv(tmp_path):
    path = tmp_path / "temp.csv"
    col = [1, 2, 3, 4, 5]
    df = pd.DataFrame({'col': col})
    df.to_csv(path, index=False)
    return str(path)
