import deepchem as dc
import numpy as np
import pandas as pd
import pytest

from deepchem_server.core import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.primitives.partition import partition


def test_partition_disk_dataset(disk_datastore):
    config.set_datastore(disk_datastore)
    dataset = dc.data.NumpyDataset(X=np.random.rand(10, 5), y=np.random.rand(10, 1))
    card = DataCard(address='', file_type='dir', data_type='DiskDataset')
    dataset_address = disk_datastore.upload_data_from_memory(dataset, "disk_data", card)

    partition_addresses = partition(dataset_address=dataset_address, n_partition=4, shuffle=False)

    assert len(partition_addresses) == 4
    total_rows = 0
    for address in partition_addresses:
        part_dataset = disk_datastore.get(address)
        total_rows += part_dataset.get_shape()[0][0]
    assert total_rows == 10

    updated_parent_card = disk_datastore.get_card(dataset_address, kind='data')
    assert updated_parent_card.n_partition == 4


def test_partition_csv_dataframe(disk_datastore):
    config.set_datastore(disk_datastore)
    df = pd.DataFrame({"smiles": [f"C{i}" for i in range(10)], "y": list(range(10))})
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    dataset_address = disk_datastore.upload_data_from_memory(df, "raw.csv", card)

    partition_addresses = partition(dataset_address=dataset_address, n_partition=4)

    assert len(partition_addresses) == 4
    total_rows = 0
    for address in partition_addresses:
        part_df = disk_datastore.get(address)
        total_rows += part_df.shape[0]
    assert total_rows == 10


def test_partition_csv_shuffle_not_supported(disk_datastore):
    config.set_datastore(disk_datastore)
    df = pd.DataFrame({"smiles": [f"C{i}" for i in range(5)], "y": list(range(5))})
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    dataset_address = disk_datastore.upload_data_from_memory(df, "raw2.csv", card)

    with pytest.raises(ValueError, match="does not support shuffling"):
        partition(dataset_address=dataset_address, n_partition=2, shuffle=True)
