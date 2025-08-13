import ast
import json
import os
import tempfile
from typing import List

import deepchem as dc
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve

from deepchem_server.core import config
from deepchem_server.core.address import DeepchemAddress
from deepchem_server.core.cards import DataCard

deepchem_server_metrics = {
    'pearson_r2_score':
        dc.metrics.Metric(dc.metrics.pearson_r2_score),
    'jaccard_score':
        dc.metrics.Metric(dc.metrics.jaccard_score),
    'prc_auc_score':
        dc.metrics.Metric(dc.metrics.prc_auc_score),
    'roc_auc_score':
        dc.metrics.Metric(dc.metrics.roc_auc_score),
    'rms_score':
        dc.metrics.Metric(dc.metrics.rms_score),
    'mae_error':
        dc.metrics.Metric(dc.metrics.mae_score),
    'bedroc_score':
        dc.metrics.Metric(dc.metrics.bedroc_score),
    'accuracy_score':
        dc.metrics.Metric(dc.metrics.accuracy_score),
    'balanced_accuracy_score':
        dc.metrics.Metric(dc.metrics.balanced_accuracy_score),
}


def prc_auc_curve(y_true: np.ndarray, y_preds: np.ndarray) -> pd.DataFrame:
    """
    Generate precision recall dataframe

    Parameters
    ----------
    y_true: np.ndarray
        true values from the dataset
    y_preds: np.ndarray
        model predictions based on the dataset

    Return
    ------
    df: pd.DataFrame
        precision recall dataframe
    """
    precision, recall, thresholds = precision_recall_curve(y_true,
                                                           y_preds[:, 1],
                                                           pos_label=1)
    thresholds = np.append(thresholds, None)  # type: ignore
    df = pd.DataFrame({
        'precision': precision,
        'recall': recall,
        'thresholds': thresholds,
    })
    return df


def model_evaluator(dataset_addresses: List[str],
                    model_address: str,
                    metrics: List[str],
                    output_key: str,
                    is_metric_plots: bool = False):
    """
    Evaluate models using featurized datasets

    Parameters
    ----------
    dataset_addresses: List[str]
        List of dataset to evaluate the model
    model_address: str
        Deepchem address of the model
    metrics: List[str]
        List of metrics to evaluate the model
    output_key: str
        Name of the evaluated output
    is_metric_plots: bool
        Whether plot based metric is used or not
    """
    # FIXME ast.literal_eval might cause security risk since it allows
    # execution of arbitary commands. We have to add literal eval parser
    # in Lark parser.
    if isinstance(dataset_addresses, str):
        dataset_addresses = ast.literal_eval(dataset_addresses)
    if isinstance(metrics, str):
        metrics = ast.literal_eval(metrics)
    if isinstance(is_metric_plots, str):
        is_metric_plots = ast.literal_eval(is_metric_plots)
    plot_metrics = {"prc_auc_curve": prc_auc_curve}
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")
    model = datastore.get(model_address, kind='model')

    if not is_metric_plots:
        eval_metrics = [
            deepchem_server_metrics[metric]
            for metric in metrics
            if metric not in plot_metrics.keys()
        ]
        if len(eval_metrics) == 0:
            raise ValueError("No non-plot metric provided to evaluate")
        dataset_scores = dict()
        for dataset_address in dataset_addresses:
            dataset = datastore.get(dataset_address)
            scores = model.evaluate(dataset, eval_metrics)
            dataset_scores[dataset_address] = scores
        dataset_scores_final: str = json.dumps(dataset_scores)
        description = "evaluation of %s datasets using %s metrics" % (
            ' '.join(dataset_addresses), ' '.join(metrics))
        card = DataCard(address='',
                        file_type='json',
                        data_type='json',
                        description=description)
        if not output_key.endswith('.json'):
            output_key = output_key + '.json'
        output_address = datastore.upload_data_from_memory(
            dataset_scores_final, DeepchemAddress.get_key(output_key), card)
    else:
        if len(metrics) > 1:
            raise Exception("only one plot metric supported")
        if metrics[0] not in plot_metrics.keys():
            raise ValueError("No plot metric provided to evaluate")
        if len(dataset_addresses) > 1:
            raise Exception("only one dataset supported for plot metric")
        dataset = datastore.get(dataset_addresses[0])
        plt_metric = plot_metrics[metrics[0]]
        y_true = dataset.y
        y_preds = model.predict(dataset)
        df = plt_metric(y_true, y_preds)

        if not output_key.endswith('.csv'):
            output_key = output_key + '.csv'

        tempdir = tempfile.TemporaryDirectory()
        temp_output_path = os.path.join(tempdir.name, output_key)
        df.to_csv(temp_output_path, index=False)

        card = DataCard(address='', file_type='csv', data_type='DataFrame')
        output_address = datastore.upload_data(
            DeepchemAddress.get_key(output_key), temp_output_path, card)
    return output_address
