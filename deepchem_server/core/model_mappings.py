import logging
import deepchem as dc
from functools import wraps
from datetime import datetime
from typing import Dict, Optional
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from deepchem_server.core.model_config_mapper import DeepChemModelConfigMapper, ModelAddressWrapper


logger = logging.getLogger(__name__)


def sklearn_model(model):
    # The wrapper here is used for distinguishing SklearnModel parameters and scikit-learn model parameters
    @wraps(model)
    def initialize_sklearn_model(model_dir: Optional[str] = None, **kwargs):
        if model_dir is None:
            return dc.models.SklearnModel(model(**kwargs))
        else:
            return dc.models.SklearnModel(model(**kwargs), model_dir=model_dir)

    return initialize_sklearn_model


model_address_map = ModelAddressWrapper({
    'random_forest_classifier':
        DeepChemModelConfigMapper(
            model_class=sklearn_model(RandomForestClassifier),
            required_init_params=None,
            optional_init_params=[
                "n_estimators", "criterion", "max_depth", "min_samples_split",
                "min_samples_leaf", "min_weight_fraction_leaf", "max_features",
                "max_leaf_nodes", "min_impurity_decrease", "bootstrap",
                "oob_score", "n_jobs", "random_state", "verbose", "warm_start",
                "class_weight", "ccp_alpha", "max_samples"
            ],
            required_train_params=None,
            optional_train_params=["sample_weight"]),
    'random_forest_regressor':
        DeepChemModelConfigMapper(
            model_class=sklearn_model(RandomForestRegressor),
            required_init_params=None,
            optional_init_params=[
                "n_estimators", "criterion", "max_depth", "min_samples_split",
                "min_samples_leaf", "min_weight_fraction_leaf", "max_features",
                "max_leaf_nodes", "min_impurity_decrease", "bootstrap",
                "oob_score", "n_jobs", "random_state", "verbose", "warm_start",
                "ccp_alpha", "max_samples"
            ],
            required_train_params=None,
            optional_train_params=["sample_weight"]),
})

LOGS = {}

def update_logs(log_error: ImportError):
    """
    Function to update logs during import errors

    Parameters
    ----------
    log_error: ImportError
        Import error object

    Example
    -------
    >>> from deepchem_server.core import model_mappings
    >>> model_mappings.LOGS == {}
    True
    >>> e = ImportError('cannot import DummyModule')
    >>> model_mappings.update_logs(e)
    >>> list(model_mappings.LOGS.values())[0]
    ImportError('cannot import DummyModule')
    """
    current_date_time = str(datetime.now())
    LOGS[current_date_time] = log_error

MODEL_FEAT_MAP = {}

try:
    # torch models
    from deepchem.models import GCNModel
    model_address_map["gcn"] = DeepChemModelConfigMapper(
        model_class=GCNModel,
        required_init_params=["n_tasks"],
        optional_init_params=[
            "graph_conv_layers", "activation", "residual", "batchnorm",
            "dropout", "predictor_hidden_feats", "predictor_dropout", "mode",
            "number_atom_features", "n_classes", "self_loop", "output_types",
            "batch_size", "learning_rate", "optimizer", "tensorboard", "wandb",
            "log_frequency", "device", "regularization_loss", "wandb_logger"
        ],
        required_train_params=None,
        optional_train_params=[
            "nb_epoch", "max_checkpoints_to_keep", "checkpoint_interval",
            "deterministic", "restore", "variables", "loss", "callbacks",
            "all_losses"
        ])
    MODEL_FEAT_MAP["gcn"] = "molgraphconv"
except ImportError as e2:
    update_logs(e2)
    logger.error(f"torch models not imported: {e2}")

model_param_map = {}
"""Creates a parameter map for models using DeepChemModelConfigMapper.

For each valid model, it retrieves initialization and training parameters
(both required and optional) and stores them in `model_param_map`. The resulting
`model_param_map` will be used to return these parameters as part of the API
endpoint's response.
"""
model_names = model_address_map.keys()

for model_name in model_names:
    if model_name == 'graphconv' or model_name == 'weave' or model_name == 'pagtn' or model_name == 'ensemble':
        continue
    model_class = str(DeepChemModelConfigMapper(
        (model_address_map[model_name])))
    model_config_map = model_address_map.get_model_config(model_class,
                                                          kind="class_name")
    if model_config_map is not None:
        model_param_map[model_name] = {
            'required_init_params':
                model_config_map.get_init_params('required'),
            'optional_init_params':
                model_config_map.get_init_params('optional'),
            'required_train_params':
                model_config_map.get_train_params('required'),
            'optional_train_params':
                model_config_map.get_train_params('optional')
        }
