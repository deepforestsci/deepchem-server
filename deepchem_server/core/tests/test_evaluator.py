import deepchem as dc
import numpy as np
import pandas as pd

from deepchem_server.core import cards, config, evaluator
from deepchem_server.core.train import train


def test_model_evaluator(disk_datastore):
    config.set_datastore(disk_datastore)
    dataset = dc.data.NumpyDataset(X=np.random.rand(100, 5),
                                   y=np.random.rand(100, 1))

    card = cards.DataCard(address='',
                          file_type='dir',
                          data_type='DiskDataset',
                          description='a disk dataset')
    dataset_address = disk_datastore.upload_data_from_memory(
        dataset, 'test_dataset', card)

    model_address = train(model_type='linear_regression',
                          dataset_address=dataset_address,
                          model_name='test_model')

    output_address = evaluator.model_evaluator(
        dataset_addresses=[dataset_address],
        model_address=model_address,
        metrics=['pearson_r2_score'],
        output_key='model_eval')

    eval_result = disk_datastore.get(output_address)
    assert dataset_address in eval_result.keys()
    assert eval_result[dataset_address]['pearson_r2_score'] > 0


def test_model_evaluator_nested_full_address(disk_datastore):
    config.set_datastore(disk_datastore)
    dataset = dc.data.NumpyDataset(X=np.random.rand(100, 5),
                                   y=np.random.rand(100, 1))

    card = cards.DataCard(address='',
                          file_type='dir',
                          data_type='DiskDataset',
                          description='a disk dataset')
    dataset_address = disk_datastore.upload_data_from_memory(
        dataset, 'test eval/test_dataset', card)
    assert dataset_address == 'deepchem://test/user/test eval/test_dataset'

    model_address = train(
        model_type='linear_regression',
        dataset_address=dataset_address,
        model_name='deepchem://test/user/test eval/test_model')
    assert model_address == 'deepchem://test/user/test eval/test_model'

    output_address = evaluator.model_evaluator(
        dataset_addresses=[dataset_address],
        model_address=model_address,
        metrics=['pearson_r2_score'],
        output_key='deepchem://test/user/test eval/model_eval')
    assert output_address == 'deepchem://test/user/test eval/model_eval.json'

    eval_result = disk_datastore.get(output_address)
    assert dataset_address in eval_result.keys()
    assert eval_result[dataset_address]['pearson_r2_score'] > 0


def test_model_evaluator_prc_auc(disk_datastore):
    """
    test model evaluator PRC AUC curve
    """
    config.set_datastore(disk_datastore)
    dataset = dc.data.NumpyDataset(X=np.random.rand(100, 5),
                                   y=np.random.randint(2, size=100))

    card = cards.DataCard(address='',
                          file_type='dir',
                          data_type='DiskDataset',
                          description='a disk dataset')
    dataset_address = disk_datastore.upload_data_from_memory(
        dataset, 'test_dataset', card)

    model_address = train(model_type='random_forest_classifier',
                          dataset_address=dataset_address,
                          model_name='test_model')

    output_address = evaluator.model_evaluator(
        dataset_addresses=[dataset_address],
        model_address=model_address,
        metrics=['prc_auc_curve'],
        output_key='model_eval_prc',
        is_metric_plots=True)

    prc_df = disk_datastore.get(output_address)
    assert isinstance(prc_df, pd.DataFrame)
    assert list(prc_df.columns) == ['precision', 'recall', 'thresholds']
    assert not (prc_df['precision'] > 1).any()
    assert not (prc_df['recall'] > 1).any()
    assert not (prc_df['thresholds'] > 1).any()
    assert not (prc_df['precision'].isna()).any()
    assert not (prc_df['recall'].isna()).any()


def test_model_evaluator_error_check(disk_datastore):
    """
    test model evaluator PRC AUC curve
    """
    config.set_datastore(disk_datastore)
    dataset = dc.data.NumpyDataset(X=np.random.rand(100, 5),
                                   y=np.random.randint(2, size=100))

    card = cards.DataCard(address='',
                          file_type='dir',
                          data_type='DiskDataset',
                          description='a disk dataset')
    dataset_address = disk_datastore.upload_data_from_memory(
        dataset, 'test_dataset', card)

    model_address = train(model_type='random_forest_classifier',
                          dataset_address=dataset_address,
                          model_name='test_model')

    try:
        evaluator.model_evaluator(dataset_addresses=[dataset_address],
                                  model_address=model_address,
                                  metrics=['prc_auc_curve'],
                                  output_key='model_eval_prc',
                                  is_metric_plots=False)
        raise Exception("Evaluator exception on non plot metric failed")
    except ValueError as e:
        assert e.__str__() == "No non-plot metric provided to evaluate"

    try:
        evaluator.model_evaluator(dataset_addresses=[dataset_address],
                                  model_address=model_address,
                                  metrics=['pearson_r2_score'],
                                  output_key='model_eval_prc',
                                  is_metric_plots=True)
        raise Exception("Evaluator exception on plot metric failed")
    except ValueError as e:
        assert e.__str__() == "No plot metric provided to evaluate"
