import pandas as pd
import pytest

from deepchem_server.core.datastore import DiskDataStore


@pytest.fixture
def disk_datastore(tmpdir):
    fds = DiskDataStore(profile_name='test',
                        project_name='user',
                        basedir=tmpdir)
    return fds


@pytest.fixture
def alternate_disk_datastore(tmpdir):
    fds = DiskDataStore(profile_name='alternate-test',
                        project_name='alternate-user',
                        basedir=tmpdir)
    return fds


@pytest.fixture
def tmpfile(tmpdir):
    text_file = tmpdir.join("temp.txt")
    text_file.write("hello world")
    return text_file


@pytest.fixture
def tmp_csv(tmpdir):
    path = tmpdir.join("temp.csv")
    col = [1, 2, 3, 4, 5]
    df = pd.DataFrame({'col': col})
    df.to_csv(path, index=False)
    return str(path)
