import pandas as pd
import pytest

from deepchem_server.utils import _init_datastore


@pytest.fixture
def disk_datastore():
    datastore = _init_datastore(profile_name="test", project_name="user", backend="local")
    return datastore


@pytest.fixture
def alternate_disk_datastore():
    datastore = _init_datastore(profile_name="alternate-test", project_name="alternate-user", backend="local")
    return datastore


@pytest.fixture
def tmp_csv(tmp_path):
    path = tmp_path / "temp.csv"
    col = [1, 2, 3, 4, 5]
    df = pd.DataFrame({'col': col})
    df.to_csv(path, index=False)
    return str(path)
