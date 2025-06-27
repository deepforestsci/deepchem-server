import os
from deepchem_server.core import config
from deepchem_server.core.feat import featurize
from deepchem_server.core.cards import DataCard
from deepchem_server.core.train import train
from deepchem_server.core.inference import infer


def test_inference(disk_datastore):
    """Test basic model inference functionality."""
    config.set_datastore(disk_datastore)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir, "assets/smiles_reg.csv")
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')

    address = disk_datastore.upload_data('example.csv', input_file, card)
    # perform ecfp featurization
    feat_address = featurize(address,
                             featurizer='ecfp',
                             output='featurized_data',
                             dataset_column='smiles',
                             label_column='log-solubility')

    # train a regression model
    model_address = train(model_type='random_forest_regressor',
                          dataset_address=feat_address,
                          model_name='ecfp_reg')

    infer_address = infer(model_address, feat_address, output='infer.csv')
    infer_df = disk_datastore.get(infer_address)
    infer_df.to_csv("test_1.csv", index=False)
    assert infer_df.shape == (10, 2)


def test_inference_nested_full_address(disk_datastore):
    """Test basic model inference functionality with nested full address."""
    config.set_datastore(disk_datastore)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir, "assets/smiles_reg.csv")
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')

    address = disk_datastore.upload_data('test infer/example.csv', input_file,
                                         card)
    assert address == 'deepchem://test/user/test infer/example.csv'

    # perform ecfp featurization
    feat_address = featurize(
        address,
        featurizer='ecfp',
        output='deepchem://test/user/test infer/featurized_data',
        dataset_column='smiles',
        label_column='log-solubility')
    assert feat_address == 'deepchem://test/user/test infer/featurized_data'

    # train a regression model
    model_address = train(model_type='random_forest_regressor',
                          dataset_address=feat_address,
                          model_name='deepchem://test/user/test infer/ecfp_reg')
    assert model_address == 'deepchem://test/user/test infer/ecfp_reg'

    infer_address = infer(model_address,
                          feat_address,
                          output='deepchem://test/user/test infer/infer.csv')
    assert infer_address == 'deepchem://test/user/test infer/infer.csv'

    infer_df = disk_datastore.get(infer_address)
    infer_df.to_csv("test_1.csv", index=False)
    assert infer_df.shape == (10, 2)


def test_featurize_and_inference(disk_datastore):
    """Test basic model inference functionality with csv files"""
    config.set_datastore(disk_datastore)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir, "assets/smiles_reg.csv")
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')

    address = disk_datastore.upload_data('example.csv', input_file, card)
    feat_address = featurize(address,
                             featurizer='molgraphconv',
                             output='gcn_feat',
                             dataset_column='smiles',
                             label_column='log-solubility')

    model_address = train(model_type='gcn',
                          dataset_address=feat_address,
                          model_name='gcn_reg',
                          init_kwargs="{'n_tasks': 1, 'mode': 'regression'}",
                          train_kwargs="{'nb_epoch': 1}")

    infer_address = infer(model_address,
                          address,
                          output='infer.csv',
                          dataset_column='smiles')
    infer_df = disk_datastore.get(infer_address)
    infer_df.to_csv("test_2.csv", index=False)
    assert infer_df.shape == (10, 2)


def test_single_datapoint_inference_rf_regress(disk_datastore):
    """Test inference for single datapoint on regression model"""
    config.set_datastore(disk_datastore)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir, "assets/smiles_reg1.csv")
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')

    address = disk_datastore.upload_data('example.csv', input_file, card)
    feat_address = featurize(address,
                             featurizer='ecfp',
                             output='rf_feat',
                             dataset_column='smiles',
                             label_column='log-solubility')

    train_input_file = os.path.join(current_dir, "assets/smiles_reg.csv")
    card2 = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')

    train_address = disk_datastore.upload_data('example.csv', train_input_file,
                                               card2)
    train_feat_address = featurize(train_address,
                                   featurizer='ecfp',
                                   output='rf_feat2',
                                   dataset_column='smiles',
                                   label_column='log-solubility')

    model_address = train(model_type='random_forest_regressor',
                          dataset_address=train_feat_address,
                          model_name='rf_reg')

    infer_address = infer(model_address,
                          feat_address,
                          output='infer.csv',
                          dataset_column='smiles')
    infer_df = disk_datastore.get(infer_address)
    infer_df.to_csv("test_2.csv", index=False)
    assert infer_df.shape == (1, 2)


def test_single_datapoint_inference_rf_class(disk_datastore):
    """Test inference for single datapoint on classification model"""
    config.set_datastore(disk_datastore)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir, "assets/sample_classifiy1.csv")
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')

    address = disk_datastore.upload_data('example.csv', input_file, card)
    feat_address = featurize(address,
                             featurizer='ecfp',
                             output='rf_feat',
                             dataset_column='smiles',
                             label_column='label')

    train_input_file = os.path.join(current_dir, "assets/sample_classifiy.csv")
    card2 = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')

    train_address = disk_datastore.upload_data('example.csv', train_input_file,
                                               card2)
    train_feat_address = featurize(train_address,
                                   featurizer='ecfp',
                                   output='rf_feat2',
                                   dataset_column='smiles',
                                   label_column='label')

    model_address = train(model_type='random_forest_classifier',
                          dataset_address=train_feat_address,
                          model_name='rf_class')

    infer_address = infer(model_address,
                          feat_address,
                          output='infer2.csv',
                          dataset_column='smiles')
    infer_df = disk_datastore.get(infer_address)
    assert infer_df.shape == (1, 3)


def test_single_datapoint_inference_gcn_class(disk_datastore):
    """Test inference for single datapoint on gcn classification model"""
    config.set_datastore(disk_datastore)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir, "assets/sample_classifiy1.csv")
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')

    address = disk_datastore.upload_data('example.csv', input_file, card)
    feat_address = featurize(address,
                             featurizer='molgraphconv',
                             output='gcn_feat',
                             dataset_column='smiles',
                             label_column='label')

    train_input_file = os.path.join(current_dir, "assets/sample_classifiy.csv")
    card2 = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')

    train_address = disk_datastore.upload_data('example.csv', train_input_file,
                                               card2)
    train_feat_address = featurize(train_address,
                                   featurizer='molgraphconv',
                                   output='gcn_feat2',
                                   dataset_column='smiles',
                                   label_column='label')

    model_address = train(
        model_type='gcn',
        dataset_address=train_feat_address,
        model_name='gcn_class',
        init_kwargs="{'n_tasks': 1, 'mode': 'classification'}",
        train_kwargs="{'nb_epoch': 1}")

    infer_address = infer(model_address,
                          feat_address,
                          output='infer2.csv',
                          dataset_column='smiles')
    infer_df = disk_datastore.get(infer_address)
    assert infer_df.shape == (1, 3)
