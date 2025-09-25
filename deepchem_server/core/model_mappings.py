from datetime import datetime
from functools import wraps
import logging
from typing import Any, Callable, Optional

import deepchem as dc
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import LinearRegression

from deepchem_server.core.model_config_mapper import DeepChemModelConfigMapper, ModelAddressWrapper


logger = logging.getLogger(__name__)


def sklearn_model(model: Callable) -> Callable:
    """Wrapper for sklearn models to integrate with DeepChem SklearnModel.

    Parameters
    ----------
    model : Callable
        A sklearn model class to be wrapped.

    Returns
    -------
    Callable
        A function that initializes a DeepChem SklearnModel with the given
        sklearn model.
    """

    # The wrapper here is used for distinguishing SklearnModel parameters and scikit-learn model parameters
    @wraps(model)
    def initialize_sklearn_model(model_dir: Optional[str] = None, **kwargs) -> Any:
        """Initialize sklearn model wrapped in DeepChem SklearnModel.

        Parameters
        ----------
        model_dir : str, optional
            Directory path for the model.
        **kwargs
            Additional keyword arguments for the sklearn model.

        Returns
        -------
        Any
            DeepChem sklearn model instance.
        """
        if model_dir is None:
            return dc.models.SklearnModel(model(**kwargs))
        else:
            return dc.models.SklearnModel(model(**kwargs), model_dir=model_dir)

    return initialize_sklearn_model


model_address_map = ModelAddressWrapper({
    "linear_regression":
        DeepChemModelConfigMapper(
            model_class=sklearn_model(LinearRegression),
            required_init_params=None,
            optional_init_params=["fit_intercept", "copy_X", "n_jobs", "positive"],
            required_train_params=None,
            optional_train_params=None,
        ),
    "random_forest_classifier":
        DeepChemModelConfigMapper(
            model_class=sklearn_model(RandomForestClassifier),
            required_init_params=None,
            optional_init_params=[
                "n_estimators",
                "criterion",
                "max_depth",
                "min_samples_split",
                "min_samples_leaf",
                "min_weight_fraction_leaf",
                "max_features",
                "max_leaf_nodes",
                "min_impurity_decrease",
                "bootstrap",
                "oob_score",
                "n_jobs",
                "random_state",
                "verbose",
                "warm_start",
                "class_weight",
                "ccp_alpha",
                "max_samples",
            ],
            required_train_params=None,
            optional_train_params=["sample_weight"],
        ),
    "random_forest_regressor":
        DeepChemModelConfigMapper(
            model_class=sklearn_model(RandomForestRegressor),
            required_init_params=None,
            optional_init_params=[
                "n_estimators",
                "criterion",
                "max_depth",
                "min_samples_split",
                "min_samples_leaf",
                "min_weight_fraction_leaf",
                "max_features",
                "max_leaf_nodes",
                "min_impurity_decrease",
                "bootstrap",
                "oob_score",
                "n_jobs",
                "random_state",
                "verbose",
                "warm_start",
                "ccp_alpha",
                "max_samples",
            ],
            required_train_params=None,
            optional_train_params=["sample_weight"],
        ),
})

LOGS = {}


def update_logs(log_error: ImportError) -> None:
    """Update logs during import errors.

    Parameters
    ----------
    log_error : ImportError
        Import error object to be logged.

    Returns
    -------
    None

    Examples
    --------
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
            "graph_conv_layers",
            "activation",
            "residual",
            "batchnorm",
            "dropout",
            "predictor_hidden_feats",
            "predictor_dropout",
            "mode",
            "number_atom_features",
            "n_classes",
            "self_loop",
            "output_types",
            "batch_size",
            "learning_rate",
            "optimizer",
            "tensorboard",
            "wandb",
            "log_frequency",
            "device",
            "regularization_loss",
            "wandb_logger",
        ],
        required_train_params=None,
        optional_train_params=[
            "nb_epoch",
            "max_checkpoints_to_keep",
            "checkpoint_interval",
            "deterministic",
            "restore",
            "variables",
            "loss",
            "callbacks",
            "all_losses",
        ],
    )
    MODEL_FEAT_MAP["gcn"] = "molgraphconv"
except ImportError as e2:
    update_logs(e2)
    logger.error(f"torch models not imported: {e2}")

model_names = model_address_map.keys()
