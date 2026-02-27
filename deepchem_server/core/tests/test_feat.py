import os
import time

import deepchem as dc
import pandas as pd
import pytest
from deepchem_server.core import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.primitives.feat import featurize


@pytest.fixture()
def zinc5k_dataset(disk_datastore):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    dataset_path = os.path.join(current_dir, "assets/zinc5k.csv")
    card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    try:
        data_address = disk_datastore.upload_data("zinc5k.csv", dataset_path, card)
    except FileExistsError:
        data_address = disk_datastore.get("zinc5k.csv")
    return data_address


def test_featurize(disk_datastore):
    """Test basic featurization functionality."""
    df = pd.DataFrame(["CCC", "CCCCC"], columns=["smiles"])
    # TODO: We should find a way to use default login for this.
    config.set_datastore(disk_datastore)
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    data_address = disk_datastore.upload_data_from_memory(df, "test.csv", card)
    dataset_address = featurize(data_address, featurizer="ecfp", output="feat_test", dataset_column="smiles")
    data = disk_datastore.get(dataset_address)
    assert data.X.shape == (2, 2048)
    card = disk_datastore.get(dataset_address + '.cdc')
    assert card.featurizer == 'ecfp'


def test_featurize_nested_full_address(disk_datastore):
    """Test basic featurization functionality."""
    df = pd.DataFrame(["CCC", "CCCCC"], columns=["smiles"])
    config.set_datastore(disk_datastore)
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    data_address = disk_datastore.upload_data_from_memory(df, "test feat/test.csv", card)
    assert data_address == "deepchem://test/user/test feat/test.csv"

    dataset_address = featurize(data_address,
                                featurizer="ecfp",
                                output="deepchem://test/user/test feat/feat_test",
                                dataset_column="smiles")
    assert dataset_address == "deepchem://test/user/test feat/feat_test"

    data = disk_datastore.get(dataset_address)
    assert data.X.shape == (2, 2048)
    card = disk_datastore.get(dataset_address + '.cdc')
    assert card.featurizer == 'ecfp'


def test_featurize_multicore(disk_datastore, zinc5k_dataset):
    """
    Test basic multicore featurization functionality with csv file.
    """
    config.set_datastore(disk_datastore)
    dataset_address = featurize(
        zinc5k_dataset,
        featurizer="ecfp",
        output="feat_test",
        dataset_column="smiles",
        label_column="logp",
        single_core_threshold=0.01,
        n_core=4,
    )  # set small `single_core_threshold` to trigger multicore feature
    data = disk_datastore.get(dataset_address)
    assert data.X.shape == (5000, 2048)
    card = disk_datastore.get(dataset_address + '.cdc')
    assert card.featurizer == 'ecfp'


def test_featurize_multicore_restart(disk_datastore, zinc5k_dataset):
    """
    Test restart functionality of multicore featurization with csv file.
    """
    config.set_datastore(disk_datastore)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    partial_checkpoints_dir_path = os.path.join(current_dir, "assets/feat_test.partial")

    output_key = f"feat_test_{time.time()}"
    checkpoint_output_key = output_key + ".partial"
    card = DataCard(
        address="",
        file_type="dir",
        data_type="dc.data.DiskDataset",
        description="partial checkpoints",
    )
    disk_datastore.upload_data(checkpoint_output_key, partial_checkpoints_dir_path, card)

    _objects = disk_datastore.list_data(include_card_files=True)

    assert f"{checkpoint_output_key}/_checkpoint/part_0_of_3.cdc" in _objects
    assert f"{checkpoint_output_key}/_checkpoint/part_1_of_3.cdc" in _objects

    dataset_address = featurize(
        zinc5k_dataset,
        featurizer="ecfp",
        output=output_key,
        dataset_column="smiles",
        label_column="logp",
        single_core_threshold=0.01,
    )  # set small `single_core_threshold` to trigger multicore feature
    data = disk_datastore.get(dataset_address)
    assert data.X.shape == (5000, 2048)
    card = disk_datastore.get(dataset_address + '.cdc')
    assert card.featurizer == 'ecfp'

    # check if checkpoint folder is deleted after the restarted job is complete
    _objects_after = disk_datastore.list_data()
    assert checkpoint_output_key + "/" not in _objects_after


def test_featurize_multicore_sdf(disk_datastore):
    """
    Test basic multicore featurization functionality with sdf file.
    """
    config.set_datastore(disk_datastore)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    dataset_path = os.path.join(current_dir, "assets/membrane_permeability.sdf")
    card = DataCard(address='', file_type='sdf', data_type='text/plain')
    data_address = disk_datastore.upload_data("membrane.sdf", dataset_path, card)
    dataset_address = featurize(data_address,
                                featurizer="ecfp",
                                output="feat_test_sdf",
                                dataset_column="None",
                                label_column="LogP(RRCK)",
                                n_core=4,
                                single_core_threshold=0.01)
    data = disk_datastore.get(dataset_address)
    assert isinstance(data, dc.data.DiskDataset)
    card = disk_datastore.get(dataset_address + '.cdc')
    assert card.featurizer == 'ecfp'


def test_featurize_sdf(disk_datastore):
    """Test basic featurization functionality."""
    # TODO: We should find a way to use default login for this.
    config.set_datastore(disk_datastore)
    card_sdf = DataCard(address='', file_type='sdf', data_type='text/plain')

    current_dir = os.path.dirname(os.path.abspath(__file__))
    sdf_test_file = os.path.join(current_dir, "assets/membrane_permeability.sdf")

    data_address = disk_datastore.upload_data("membrane.sdf", sdf_test_file, card_sdf)
    dataset_address = featurize(data_address,
                                featurizer="ecfp",
                                output="feat_test",
                                dataset_column="None",
                                label_column="LogP(RRCK)")
    data = disk_datastore.get(dataset_address)

    assert data.X.shape == (10, 2048)
    card = disk_datastore.get(dataset_address + '.cdc')
    assert card.featurizer == 'ecfp'


def test_featurize_with_label(disk_datastore):
    """Test that featurization works with labels."""
    df = pd.DataFrame([["CCC", 0], ["CCCCC", 1]], columns=["smiles", "label"])
    config.set_datastore(disk_datastore)
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    data_address = disk_datastore.upload_data_from_memory(df, "test.csv", card)
    dataset_address = featurize(data_address,
                                featurizer="ecfp",
                                output="feat_test",
                                dataset_column="smiles",
                                label_column="label")
    data = disk_datastore.get(dataset_address)
    assert data.X.shape == (2, 2048)
    assert data.y.shape == (2, 1)

    card = disk_datastore.get(dataset_address + '.cdc')
    assert card.featurizer == 'ecfp'


def test_featurize_with_feat_kwargs(disk_datastore):
    """Test featurize method with feat_kwargs"""
    df = pd.DataFrame([["CCC", 0], ["CCCCC", 1]], columns=["smiles", "label"])
    config.set_datastore(disk_datastore)
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    data_address = disk_datastore.upload_data_from_memory(df, "test.csv", card)
    feat_kwargs = {'size': 1024}
    dataset_address = featurize(data_address,
                                featurizer="ecfp",
                                output="feat_test",
                                dataset_column="smiles",
                                label_column="label",
                                feat_kwargs=feat_kwargs)
    data = disk_datastore.get(dataset_address)
    assert data.X.shape == (2, 1024)
    assert data.y.shape == (2, 1)

    card = disk_datastore.get(dataset_address + '.cdc')
    assert card.featurizer == 'ecfp'
    assert card.feat_kwargs == {'size': 1024}


def test_molgraphconv_feat(disk_datastore):
    """
    Test molgraphconv featurizer functionality.
    """
    df = pd.DataFrame(["C1=CC=CN=C1", "O=C(NCc1cc(OC)c(O)cc1)CCCC/C=C/C(C)C"], columns=["smiles"])
    config.set_datastore(disk_datastore)
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    data_address = disk_datastore.upload_data_from_memory(df, "test.csv", card)
    dataset_address = featurize(data_address, featurizer="molgraphconv", output="feat_test", dataset_column="smiles")
    data = disk_datastore.get(dataset_address)
    assert data.X.shape == (2,)
    assert isinstance(data.X[0], dc.feat.GraphData)

    # assert "C1=CC=CN=C1"
    assert data.X[0].num_nodes == 6
    assert data.X[0].num_node_features == 30
    assert data.X[0].num_edges == 12

    card = disk_datastore.get(dataset_address + '.cdc')
    assert card.featurizer == 'molgraphconv'
