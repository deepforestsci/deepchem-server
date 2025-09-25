import deepchem as dc
import pandas as pd

from deepchem_server.core import config, featurize, train
from deepchem_server.core.cards import DataCard


def test_train(disk_datastore):
    """Test basic model training functionality."""
    config.set_datastore(disk_datastore)
    df = pd.DataFrame([["CCC", 0], ["CCCCC", 1]], columns=["smiles", "label"])
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    data_address = disk_datastore.upload_data_from_memory(df, "test.csv", card)
    dataset_address = featurize(data_address,
                                featurizer="ecfp",
                                output="feat_test",
                                dataset_column="smiles",
                                label_column="label")

    model_address = train(model_type='random_forest_classifier',
                          dataset_address=dataset_address,
                          model_name='rf_classifier',
                          init_kwargs={},
                          train_kwargs={})

    model = config.get_datastore().get_model(model_address)
    assert isinstance(model, dc.models.Model)


def test_train_nested_full_address(disk_datastore):
    """Test basic model training functionality."""
    config.set_datastore(disk_datastore)
    df = pd.DataFrame([["CCC", 0], ["CCCCC", 1]], columns=["smiles", "label"])
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    data_address = disk_datastore.upload_data_from_memory(df, "test train/test.csv", card)
    assert data_address == "deepchem://test/user/test train/test.csv"

    dataset_address = featurize(data_address,
                                featurizer="ecfp",
                                output="deepchem://test/user/test train/feat_test",
                                dataset_column="smiles",
                                label_column="label")
    assert dataset_address == "deepchem://test/user/test train/feat_test"

    model_address = train(model_type='random_forest_classifier',
                          dataset_address=dataset_address,
                          model_name='deepchem://test/user/test train/rf_classifier',
                          init_kwargs={},
                          train_kwargs={})

    model = config.get_datastore().get_model(model_address)
    assert isinstance(model, dc.models.Model)
    assert model_address == "deepchem://test/user/test train/rf_classifier"
