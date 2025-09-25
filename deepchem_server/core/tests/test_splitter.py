import os

import deepchem as dc
import numpy as np
import pandas as pd

from deepchem_server.core import cards, config, featurize, splitter
from deepchem_server.core.cards import DataCard


def test_k_fold_split_deepchem_dataset(disk_datastore):
    dataset = dc.data.NumpyDataset(X=np.random.rand(100, 5), y=np.random.rand(100, 1))

    card = cards.DataCard(address='', file_type='dir', data_type='DiskDataset', description='bace.csv')
    address = disk_datastore.upload_data_from_memory(dataset, 'dd', card)
    config.set_datastore(disk_datastore)
    addresses = splitter.k_fold_split(splitter_type='index', dataset_address=address, k=4)
    n_test_data_points = 0
    for train_address, test_address in addresses:
        test_data = disk_datastore.get(test_address)
        n_test_data_points += test_data.X.shape[0]
    assert n_test_data_points == 100


def test_k_fold_scaffold_split_deepchem_dataset(disk_datastore):
    """Test k-fold scaffold split on deepchem dataset"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    dataset_path = os.path.join(current_dir, "assets/zinc5k.csv")
    df = pd.read_csv(dataset_path)
    config.set_datastore(disk_datastore)
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    data_address = disk_datastore.upload_data_from_memory(df, "test.csv", card)
    dataset_address = featurize(data_address,
                                featurizer="ecfp",
                                output="feat",
                                dataset_column="smiles",
                                label_column="logp")
    addresses = splitter.k_fold_split(splitter_type='scaffold', dataset_address=dataset_address, k=3)
    n_test_data_points = 0
    for train_address, test_address in addresses:
        test_data = disk_datastore.get(test_address)
        n_test_data_points += test_data.X.shape[0]
        dataset_card = disk_datastore.get_card(train_address, kind='data')
        assert dataset_card.featurizer == "ecfp"
        dataset_card = disk_datastore.get_card(test_address, kind='data')
        assert dataset_card.featurizer == "ecfp"
    assert n_test_data_points == 5000


def test_k_fold_stratified_split_deepchem_dataset(disk_datastore):
    """Test k-fold stratified split on deepchem dataset"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    dataset_path = os.path.join(current_dir, "assets/zinc5k.csv")
    df = pd.read_csv(dataset_path)
    config.set_datastore(disk_datastore)
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    data_address = disk_datastore.upload_data_from_memory(df, "test.csv", card)
    dataset_address = featurize(data_address,
                                featurizer="ecfp",
                                output="feat",
                                dataset_column="smiles",
                                label_column="logp")
    addresses = splitter.k_fold_split(splitter_type='random_stratified', dataset_address=dataset_address, k=3)
    n_test_data_points = 0
    for train_address, test_address in addresses:
        test_data = disk_datastore.get(test_address)
        n_test_data_points += test_data.X.shape[0]
        dataset_card = disk_datastore.get_card(train_address, kind='data')
        assert dataset_card.featurizer == "ecfp"
        dataset_card = disk_datastore.get_card(test_address, kind='data')
        assert dataset_card.featurizer == "ecfp"
    assert n_test_data_points == 5000


def test_tvt_scaffold_split_deepchem_dataset(disk_datastore):
    """Test train-valid-test scaffold split on deepchem dataset"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    dataset_path = os.path.join(current_dir, "assets/zinc5k.csv")
    df = pd.read_csv(dataset_path)
    config.set_datastore(disk_datastore)
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    data_address = disk_datastore.upload_data_from_memory(df, "test.csv", card)
    dataset_address = featurize(data_address,
                                featurizer="ecfp",
                                output="feat",
                                dataset_column="smiles",
                                label_column="logp")
    config.set_datastore(disk_datastore)
    addresses = splitter.train_valid_test_split(
        splitter_type='scaffold',
        frac_train=0.7,
        frac_valid=0.2,
        frac_test=0.1,
        dataset_address=dataset_address,
    )
    assert len(addresses) == 3
    train_address, valid_address, test_address = addresses
    train_data = disk_datastore.get(train_address)
    valid_data = disk_datastore.get(valid_address)
    test_data = disk_datastore.get(test_address)

    train_dataset_card = disk_datastore.get_card(train_address, kind='data')
    valid_dataset_card = disk_datastore.get_card(valid_address, kind='data')
    test_dataset_card = disk_datastore.get_card(test_address, kind='data')
    assert train_dataset_card.featurizer == "ecfp"
    assert valid_dataset_card.featurizer == "ecfp"
    assert test_dataset_card.featurizer == "ecfp"

    assert len(train_data) == 3500
    assert len(valid_data) == 1000
    assert len(test_data) == 500


def test_tvt_stratified_split_deepchem_dataset(disk_datastore):
    """Test train-valid-test stratified split on deepchem dataset"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    dataset_path = os.path.join(current_dir, "assets/zinc5k.csv")
    df = pd.read_csv(dataset_path)
    config.set_datastore(disk_datastore)
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    data_address = disk_datastore.upload_data_from_memory(df, "test.csv", card)
    dataset_address = featurize(data_address,
                                featurizer="ecfp",
                                output="feat",
                                dataset_column="smiles",
                                label_column="logp")
    config.set_datastore(disk_datastore)
    addresses = splitter.train_valid_test_split(
        splitter_type='random_stratified',
        frac_train=0.7,
        frac_valid=0.2,
        frac_test=0.1,
        dataset_address=dataset_address,
    )
    assert len(addresses) == 3
    train_address, valid_address, test_address = addresses
    train_data = disk_datastore.get(train_address)
    valid_data = disk_datastore.get(valid_address)
    test_data = disk_datastore.get(test_address)

    train_dataset_card = disk_datastore.get_card(train_address, kind='data')
    valid_dataset_card = disk_datastore.get_card(valid_address, kind='data')
    test_dataset_card = disk_datastore.get_card(test_address, kind='data')
    assert train_dataset_card.featurizer == "ecfp"
    assert valid_dataset_card.featurizer == "ecfp"
    assert test_dataset_card.featurizer == "ecfp"

    assert len(train_data) == 3500
    assert len(valid_data) == 1000
    assert len(test_data) == 500
