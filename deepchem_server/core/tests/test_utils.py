import os

import deepchem as dc
import numpy as np
import pytest

from deepchem_server.core.cards import DataCard
from deepchem_server.utils import _download_data, _init_datastore


def test_download_data_file(disk_datastore, tmp_csv):
    """Test _download_data helper for a file."""
    data_card = DataCard(
        address="",
        file_type="csv",
        data_type="pandas.DataFrame",
        description="this is a pandas dataframe",
    )
    data_address = disk_datastore.upload_data("test_download.csv", tmp_csv, data_card)

    file_path, is_directory, object_name = _download_data(profile_name="test",
                                                          project_name="user",
                                                          address=data_address,
                                                          backend="local")

    assert object_name == "test_download.csv"
    assert is_directory is False
    assert os.path.isfile(file_path)


def test_download_data_directory(disk_datastore):
    """Test _download_data helper for a directory."""
    X = np.random.rand(10, 10)
    y = np.random.rand(10)
    data_card = DataCard(
        address="",
        file_type="dir",
        data_type="dc.data.DiskDataset",
        description="this is a disk dataset",
    )
    data = dc.data.NumpyDataset(X, y)
    data_address = disk_datastore.upload_data_from_memory(data, "test_download_dataset", data_card)

    file_path, is_directory, object_name = _download_data(profile_name="test",
                                                          project_name="user",
                                                          address=data_address,
                                                          backend="local")

    assert object_name == "test_download_dataset"
    assert is_directory is True
    assert os.path.isdir(file_path)


def test_download_data_not_found():
    """Test _download_data helper for non-existent object."""
    non_existent_address = "deepchem://test/user/non_existent.csv"

    with pytest.raises(FileNotFoundError) as exc_info:
        _download_data(profile_name="test", project_name="user", address=non_existent_address, backend="local")

    assert "Object not found at address" in str(exc_info.value)


def test_download_data_nested_path(disk_datastore, tmp_csv):
    """Test _download_data helper for nested path."""
    data_card = DataCard(
        address="",
        file_type="csv",
        data_type="pandas.DataFrame",
        description="this is a pandas dataframe",
    )
    data_address = disk_datastore.upload_data("nested/path/test.csv", tmp_csv, data_card)

    file_path, is_directory, object_name = _download_data(profile_name="test",
                                                          project_name="user",
                                                          address=data_address,
                                                          backend="local")

    assert object_name == "test.csv"
    assert is_directory is False
    assert os.path.isfile(file_path)


def test_init_datastore(tmp_path):
    """Test _init_datastore helper function."""
    datastore = _init_datastore(profile_name="test_profile", project_name="test_project", backend="local")

    assert datastore is not None
    assert "test_profile" in datastore.storage_loc
    assert "test_project" in datastore.storage_loc


def test_init_datastore_invalid_backend():
    """Test _init_datastore with invalid backend."""
    with pytest.raises(NotImplementedError) as exc_info:
        _init_datastore(profile_name="test", project_name="test", backend="invalid_backend")

    assert "invalid_backend backend not implemented" in str(exc_info.value)
