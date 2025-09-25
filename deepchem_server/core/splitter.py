import os
import tempfile

import deepchem as dc
import numpy as np
import pandas as pd

from deepchem_server.core import config
from deepchem_server.core.address import DeepchemAddress
from deepchem_server.core.cards import DataCard


splitter_map = {
    'random': dc.splits.RandomSplitter(),
    'index': dc.splits.IndexSplitter(),
    'scaffold': dc.splits.ScaffoldSplitter(),
    'random_stratified': dc.splits.RandomStratifiedSplitter(),
}


def _k_fold_split_dc_dataset(splitter_type, dataset_address, k):
    datastore = config.get_datastore()
    address_key = DeepchemAddress.get_key(dataset_address)
    splitter = splitter_map[splitter_type]
    dataset = datastore.get(dataset_address)
    dataset_card = datastore.get_card(dataset_address, kind='data')
    featurizer = dataset_card.featurizer
    splits = splitter.k_fold_split(dataset, k)
    addresses = []
    for i, (train, test) in enumerate(splits):
        card = DataCard(address='', file_type='dir', data_type=type(train).__name__, featurizer=featurizer)
        train_address = datastore.upload_data_from_memory(data=train,
                                                          datastore_filename=address_key + f'_train{i}',
                                                          card=card)
        test_address = datastore.upload_data_from_memory(data=test,
                                                         datastore_filename=address_key + f'_test{i}',
                                                         card=card)
        addresses.append((train_address, test_address))
    return addresses


def _k_fold_split_csv(splitter_type, dataset_address, k):
    address_key = DeepchemAddress.get_key(dataset_address).strip('.csv')
    filedir = tempfile.TemporaryDirectory()
    filepath = os.path.join(filedir.name, 'temp.csv')
    datastore = config.get_datastore()
    datastore.download_object(dataset_address, filepath)
    df = pd.read_csv(filepath)
    index = list(df.index)
    test_size = len(index) // (k + 1)
    test_indices = []
    if splitter_type == 'random':
        np.random.shuffle(index)
    elif splitter_type == 'index':
        # split using uniform indices
        pass
    for i in range(k):
        test_indices.append(index[i * test_size:(i + 1) * test_size])
    test_indices.append(index[((i + 1) * test_size):])
    addresses = []
    for i, test_index in enumerate(test_indices):
        df_train = df[~df.index.isin(test_index)]
        df_test = df[df.index.isin(test_index)]
        card = DataCard(address='', file_type='dir', data_type=type(df_train).__name__)
        train_address = datastore.upload_data_from_memory(data=df_train,
                                                          datastore_filename=address_key + f'_train{i}.csv',
                                                          card=card)
        test_address = datastore.upload_data_from_memory(data=df_test,
                                                         datastore_filename=address_key + f'_test{i}.csv',
                                                         card=card)
        addresses.append((train_address, test_address))
    return addresses


def k_fold_split(splitter_type: str, dataset_address: str, k: int):
    if isinstance(k, str):
        k = int(k)
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")
    data_card = datastore.get(dataset_address + '.cdc')
    if data_card.data_type in ['dc.data.DiskDataset', 'DiskDataset', 'dc.data.NumpyDataset', 'NumpyDataset']:
        return _k_fold_split_dc_dataset(splitter_type, dataset_address, k)
    elif data_card.data_type in ['pandas.DataFrame', 'DataFrame']:
        return _k_fold_split_csv(splitter_type, dataset_address, k)


def _train_test_valid_split_dc_dataset(splitter_type, dataset_address, frac_train, frac_valid, frac_test):
    datastore = config.get_datastore()
    address_key = DeepchemAddress.get_key(dataset_address)
    splitter = splitter_map[splitter_type]
    dataset = datastore.get(dataset_address)
    dataset_card = datastore.get_card(dataset_address, kind='data')
    featurizer = dataset_card.featurizer
    train, valid, test = splitter.train_valid_test_split(dataset,
                                                         frac_train=frac_train,
                                                         frac_valid=frac_valid,
                                                         frac_test=frac_test)
    card = DataCard(address='',
                    file_type='dir',
                    data_type=type(train).__name__,
                    description='train dataset',
                    featurizer=featurizer)
    train_address = datastore.upload_data_from_memory(data=train, datastore_filename=address_key + '_train', card=card)
    card.description = 'validation dataset'
    valid_address = datastore.upload_data_from_memory(data=valid, datastore_filename=address_key + '_valid', card=card)
    card.description = 'test dataset'
    test_address = datastore.upload_data_from_memory(data=test, datastore_filename=address_key + '_test', card=card)
    return [train_address, valid_address, test_address]


def _train_test_valid_split_csv(splitter_type, dataset_address, frac_train, frac_valid, frac_test):
    # Not suitable for splitting very large datasets
    address_key = DeepchemAddress.get_key(dataset_address).strip('.csv')
    filedir = tempfile.TemporaryDirectory()
    filepath = os.path.join(filedir.name, 'temp.csv')
    datastore = config.get_datastore()
    datastore.download_object(dataset_address, filepath)
    df = pd.read_csv(filepath)
    index = list(df.index)
    index_len = len(index)
    np.random.shuffle(index)
    train_index, valid_index, test_index = np.split(
        index, [int(frac_train * index_len), int((frac_train + frac_valid) * index_len)])
    train_df = df.iloc[train_index].reset_index(drop=True)
    test_df = df.iloc[test_index].reset_index(drop=True)
    valid_df = df.iloc[test_index].reset_index(drop=True)
    card = DataCard(address='', file_type='csv', data_type=type(train_df).__name__, description='train dataset')
    train_address = datastore.upload_data_from_memory(train_df,
                                                      datastore_filename=address_key + '_train.csv',
                                                      card=card)
    card.description = 'test dataset'
    test_address = datastore.upload_data_from_memory(test_df, datastore_filename=address_key + '_test.csv', card=card)
    card.description = 'validation dataset'
    valid_address = datastore.upload_data_from_memory(valid_df,
                                                      datastore_filename=address_key + '_valid.csv',
                                                      card=card)
    return [train_address, valid_address, test_address]


def train_valid_test_split(splitter_type: str,
                           dataset_address: str,
                           frac_train: float = 0.8,
                           frac_valid: float = 0.1,
                           frac_test: float = 0.1):
    frac_train, frac_valid, frac_test = float(frac_train), float(frac_valid), float(frac_test)
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")
    data_card = datastore.get(dataset_address + '.cdc')
    if data_card.data_type in ['dc.data.DiskDataset', 'DiskDataset', 'dc.data.NumpyDataset', 'NumpyDataset']:
        return _train_test_valid_split_dc_dataset(splitter_type, dataset_address, frac_train, frac_valid, frac_test)
    elif data_card.data_type in ['pandas.DataFrame', 'DataFrame']:
        return _train_test_valid_split_csv(splitter_type, dataset_address, frac_train, frac_valid, frac_test)
