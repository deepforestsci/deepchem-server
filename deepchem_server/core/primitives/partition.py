"""Primitive to partition datasets"""

import os
import shutil
import tempfile
from typing import List, Optional

import numpy as np
import pandas as pd

from deepchem_server.core.common import config
from deepchem_server.core.common.address import DeepchemAddress
from deepchem_server.core.common.cards import DataCard


def _partition_disk_dataset(dataset_address: str, shuffle: bool, n_partition: int) -> List[str]:
    """Partition a DeepChem DiskDataset into N partitions.

    Parameters
    ----------
    dataset_address: str
        datastore address of dataset to partition
    shuffle: bool
        Whether to shuffle before partitioning
    n_partition: int
        Number of partitions to create

    Returns
    -------
    List[str]
        List of datastore addresses of the partitioned datasets
    
    Examples
    --------
    >>> from deepchem_server.core.primitives import partition
    >>> dataset_address = "deepchem://deepchem/data/zinc5k"
    >>> partitioned_dataset_addresses = partition._partition_disk_dataset(dataset_address, shuffle=False, n_partition=4)
    >>> print(partitioned_dataset_addresses)
    ['deepchem://deepchem/data/zinc5k_partition0', 'deepchem://deepchem/data/zinc5k_partition1', 'deepchem://deepchem/data/zinc5k_partition2', 'deepchem://deepchem/data/zinc5k_partition3']
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")

    dataset = datastore.get(dataset_address)
    if dataset is None:
        raise ValueError(f"Dataset not found at address: {dataset_address}")

    key = DeepchemAddress.get_key(dataset_address)
    n_elements = dataset.get_shape()[0][0]
    indices = np.arange(n_elements)
    if shuffle:
        np.random.shuffle(indices)

    elements_per_partition = n_elements / n_partition
    partitioned_dataset_addresses: List[str] = []
    dataset_card = datastore.get_card(dataset_address, kind="data")
    featurizer = dataset_card.featurizer if dataset_card is not None else None

    with tempfile.TemporaryDirectory() as tempdir:
        for i in range(n_partition):
            partition_indices = indices[int(i * elements_per_partition):int((i + 1) * elements_per_partition)]
            partition_path = os.path.join(tempdir, f"partition{i}")
            partition_dataset = dataset.select(partition_indices, select_dir=partition_path)
            card = DataCard(
                address="",
                file_type="dir",
                data_type=type(partition_dataset).__name__,
                n_partition=n_partition,
                partition_id=i,
                parent=dataset_address,
                featurizer=featurizer,
            )
            partition_key = f"{key}partition{i}"
            address = datastore.upload_data_from_memory(
                data=partition_dataset,
                datastore_filename=partition_key,
                card=card,
            )
            partitioned_dataset_addresses.append(address)
            shutil.rmtree(partition_dataset.data_dir, ignore_errors=True)

    return partitioned_dataset_addresses


def _partition_csv_dataframe(dataset_address: str, n_partition: int) -> List[str]:
    """Partition a CSV dataframe into N CSV files.

    Parameters
    ----------
    dataset_address: str
        datastore address of dataset to partition
    n_partition: int
        Number of partitions to create

    Returns
    -------
    List[str]
        List of datastore addresses of the partitioned datasets
    
    Examples
    --------
    >>> from deepchem_server.core.primitives import partition
    >>> dataset_address = "deepchem://deepchem/data/zinc5k"
    >>> partitioned_dataset_addresses = partition._partition_csv_dataframe(dataset_address, n_partition=4)
    >>> print(partitioned_dataset_addresses)
    ['deepchem://deepchem/data/zinc5k_partition0', 'deepchem://deepchem/data/zinc5k_partition1', 'deepchem://deepchem/data/zinc5k_partition2', 'deepchem://deepchem/data/zinc5k_partition3']
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")

    key = DeepchemAddress.get_key(dataset_address)
    stem = key.rsplit(".", 1)[0] if "." in key else key

    with tempfile.TemporaryDirectory() as dataset_dir:
        dataset_path = os.path.join(dataset_dir, "dataset.csv")
        datastore.download_object(dataset_address, dataset_path)

        # subtract 1 to skip the header row
        num_rows = max(0, sum(1 for _ in open(dataset_path, "r", encoding="utf-8")) - 1)
        num_rows_per_partition = int(np.ceil(num_rows / n_partition)) if num_rows > 0 else 1
        chunksize = min(100_000, num_rows_per_partition)

        partition_file_paths: List[str] = []
        with tempfile.TemporaryDirectory() as tempdir:
            partition_id = 0
            partition_file_path = os.path.join(tempdir, f"{stem}partition{partition_id}.csv")
            row_count_in_partition = 0

            for df_block in pd.read_csv(dataset_path, chunksize=chunksize):
                future_rows = row_count_in_partition + df_block.shape[0]

                if future_rows > num_rows_per_partition:
                    excess_rows = future_rows - num_rows_per_partition

                    if row_count_in_partition == 0:
                        df_block[:-excess_rows].to_csv(partition_file_path, mode="w", index=False)
                    else:
                        df_block[:-excess_rows].to_csv(partition_file_path, mode="a", header=False, index=False)

                    if partition_file_path not in partition_file_paths:
                        partition_file_paths.append(partition_file_path)

                    partition_id += 1
                    partition_file_path = os.path.join(tempdir, f"{stem}partition{partition_id}.csv")
                    df_block[-excess_rows:].to_csv(partition_file_path, mode="w", index=False)
                    row_count_in_partition = excess_rows
                else:
                    if os.path.exists(partition_file_path):
                        df_block.to_csv(partition_file_path, mode="a", header=False, index=False)
                    else:
                        df_block.to_csv(partition_file_path, mode="w", index=False)
                    row_count_in_partition = future_rows

            if (os.path.exists(partition_file_path) and partition_file_path not in partition_file_paths):
                partition_file_paths.append(partition_file_path)

            partitioned_dataset_addresses: List[str] = []
            for i, file_path in enumerate(partition_file_paths):
                card = DataCard(
                    address="",
                    file_type="csv",
                    data_type="DataFrame",
                    n_partition=n_partition,
                    partition_id=i,
                    parent=dataset_address,
                )
                address = datastore.upload_data(
                    datastore_filename=f"{stem}partition{i}.csv",
                    filename=file_path,
                    card=card,
                )
                partitioned_dataset_addresses.append(address)

    return partitioned_dataset_addresses


def partition(dataset_address: str, shuffle: bool = False, n_partition: int = 4) -> List[str]:
    """Partition a dataset into N datasets.

    Parameters
    ----------
    dataset_address: str
        datastore address of dataset to partition
    shuffle: bool
        Whether to shuffle before partitioning
    n_partition: int
        Number of partitions to create

    Returns
    -------
    List[str]
        List of datastore addresses of the partitioned datasets

    Notes
    -----
    - Supported inputs: DeepChem ``DiskDataset`` and CSV ``DataFrame``.
    - CSV partitioning does not support shuffling.
    
    Examples
    --------
    >>> from deepchem_server.core.primitives import partition
    >>> dataset_address = "deepchem://deepchem/data/zinc5k"
    >>> partitioned_dataset_addresses = partition(dataset_address, n_partition=4)
    >>> print(partitioned_dataset_addresses)
    ['deepchem://deepchem/data/zinc5k_partition0', 'deepchem://deepchem/data/zinc5k_partition1', 'deepchem://deepchem/data/zinc5k_partition2', 'deepchem://deepchem/data/zinc5k_partition3']
    
    >>> partitioned_dataset_addresses = partition(dataset_address, n_partition=4, shuffle=True)
    >>> print(partitioned_dataset_addresses)
    ['deepchem://deepchem/data/zinc5k_partition0', 'deepchem://deepchem/data/zinc5k_partition1', 'deepchem://deepchem/data/zinc5k_partition2', 'deepchem://deepchem/data/zinc5k_partition3']
    """
    if isinstance(shuffle, str):
        shuffle = shuffle.lower() == "true"
    if isinstance(n_partition, str):
        n_partition = int(n_partition)
    if n_partition <= 0:
        raise ValueError("n_partition must be a positive integer")

    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")

    data_card: Optional[DataCard] = datastore.get_card(dataset_address, kind="data")  # type: ignore
    if data_card is None:
        raise ValueError(f"Datacard not found for dataset address: {dataset_address}")

    if data_card.data_type in ("dc.data.DiskDataset", "DiskDataset"):
        partitioned_datasets = _partition_disk_dataset(dataset_address, shuffle, n_partition)
    elif data_card.data_type in ("pandas.DataFrame", "DataFrame"):
        if shuffle:
            raise ValueError("CSV/DataFrame partitioning does not support shuffling.")
        if data_card.file_type != "csv":
            raise NotImplementedError(f"DataFrame partitioning for file_type '{data_card.file_type}' is not supported.")
        partitioned_datasets = _partition_csv_dataframe(dataset_address, n_partition)
    else:
        raise NotImplementedError(f"Dataset type '{data_card.data_type}' is not supported for partitioning.")

    # Update parent card metadata with partition information.
    data_card.update_card("n_partition", n_partition)
    parent_key = DeepchemAddress.get_key(dataset_address)
    datastore.upload_data_from_memory(data=data_card, datastore_filename=f"{parent_key}.cdc", card=None)

    return partitioned_datasets
