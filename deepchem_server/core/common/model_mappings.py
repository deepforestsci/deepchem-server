"""
Model mappings for DeepChem Server.

This module provides lazy-loaded model configurations to avoid importing
heavy ML dependencies (deepchem, sklearn, torch) until actually needed.
"""
import logging
from datetime import datetime
from functools import wraps
from typing import Any, Callable, Dict, List, Optional

from deepchem_server.core.common.model_config_mapper import (
    DeepChemModelConfigMapper,
)

logger = logging.getLogger(__name__)

LOGS: Dict[str, ImportError] = {}


def update_logs(log_error: ImportError) -> None:
    """Update logs during import errors.

    Parameters
    ----------
    log_error : ImportError
        Import error object to be logged.
    """
    current_date_time = str(datetime.now())
    LOGS[current_date_time] = log_error


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

    @wraps(model)
    def initialize_sklearn_model(model_dir: Optional[str] = None, **kwargs) -> Any:
        """Initialize sklearn model wrapped in DeepChem SklearnModel."""
        # Lazy import of deepchem
        import deepchem as dc
        
        if model_dir is None:
            return dc.models.SklearnModel(model(**kwargs))
        else:
            return dc.models.SklearnModel(model(**kwargs), model_dir=model_dir)

    return initialize_sklearn_model


class LazyModelRegistry:
    """Registry that lazy-loads model configurations.
    
    Models are only loaded when accessed, avoiding importing heavy
    ML dependencies at module load time.
    """
    
    def __init__(self) -> None:
        self._loaded_models: Dict[str, DeepChemModelConfigMapper] = {}
        self._model_loaders: Dict[str, Callable[[], DeepChemModelConfigMapper]] = {}
        self._model_feat_map: Dict[str, str] = {}
        self._initialized = False
    
    def _register_loaders(self) -> None:
        """Register lazy loaders for all models."""
        if self._initialized:
            return
        
        # Sklearn models
        self._model_loaders["linear_regression"] = self._load_linear_regression
        self._model_loaders["random_forest_classifier"] = self._load_random_forest_classifier
        self._model_loaders["random_forest_regressor"] = self._load_random_forest_regressor
        
        # PyTorch/DGL models
        self._model_loaders["gcn"] = self._load_gcn
        
        self._initialized = True
    
    def _load_linear_regression(self) -> DeepChemModelConfigMapper:
        """Lazy load LinearRegression model config."""
        from sklearn.linear_model import LinearRegression
        return DeepChemModelConfigMapper(
            model_class=sklearn_model(LinearRegression),
            required_init_params=None,
            optional_init_params=["fit_intercept", "copy_X", "n_jobs", "positive"],
            required_train_params=None,
            optional_train_params=None,
        )
    
    def _load_random_forest_classifier(self) -> DeepChemModelConfigMapper:
        """Lazy load RandomForestClassifier model config."""
        from sklearn.ensemble import RandomForestClassifier
        return DeepChemModelConfigMapper(
            model_class=sklearn_model(RandomForestClassifier),
            required_init_params=None,
            optional_init_params=[
                "n_estimators", "criterion", "max_depth", "min_samples_split",
                "min_samples_leaf", "min_weight_fraction_leaf", "max_features",
                "max_leaf_nodes", "min_impurity_decrease", "bootstrap", "oob_score",
                "n_jobs", "random_state", "verbose", "warm_start", "class_weight",
                "ccp_alpha", "max_samples",
            ],
            required_train_params=None,
            optional_train_params=["sample_weight"],
        )
    
    def _load_random_forest_regressor(self) -> DeepChemModelConfigMapper:
        """Lazy load RandomForestRegressor model config."""
        from sklearn.ensemble import RandomForestRegressor
        return DeepChemModelConfigMapper(
            model_class=sklearn_model(RandomForestRegressor),
            required_init_params=None,
            optional_init_params=[
                "n_estimators", "criterion", "max_depth", "min_samples_split",
                "min_samples_leaf", "min_weight_fraction_leaf", "max_features",
                "max_leaf_nodes", "min_impurity_decrease", "bootstrap", "oob_score",
                "n_jobs", "random_state", "verbose", "warm_start",
                "ccp_alpha", "max_samples",
            ],
            required_train_params=None,
            optional_train_params=["sample_weight"],
        )
    
    def _load_gcn(self) -> DeepChemModelConfigMapper:
        """Lazy load GCN model config."""
        try:
            from deepchem.models import GCNModel
            self._model_feat_map["gcn"] = "molgraphconv"
            return DeepChemModelConfigMapper(
                model_class=GCNModel,
                required_init_params=["n_tasks"],
                optional_init_params=[
                    "graph_conv_layers", "activation", "residual", "batchnorm",
                    "dropout", "predictor_hidden_feats", "predictor_dropout",
                    "mode", "number_atom_features", "n_classes", "self_loop",
                    "output_types", "batch_size", "learning_rate", "optimizer",
                    "tensorboard", "wandb", "log_frequency", "device",
                    "regularization_loss", "wandb_logger",
                ],
                required_train_params=None,
                optional_train_params=[
                    "nb_epoch", "max_checkpoints_to_keep", "checkpoint_interval",
                    "deterministic", "restore", "variables", "loss", "callbacks",
                    "all_losses",
                ],
            )
        except ImportError as e:
            update_logs(e)
            logger.error(f"GCN model not available: {e}")
            raise
    
    def __getitem__(self, name: str) -> Any:
        """Get a model config, loading it lazily if needed."""
        self._register_loaders()
        
        if name not in self._loaded_models:
            if name not in self._model_loaders:
                raise KeyError(f"Unknown model: {name}")
            self._loaded_models[name] = self._model_loaders[name]()
        
        return self._loaded_models[name].get_model_class()
    
    def __contains__(self, name: str) -> bool:
        """Check if a model is registered."""
        self._register_loaders()
        return name in self._model_loaders
    
    def keys(self) -> List[str]:
        """Get all registered model names."""
        self._register_loaders()
        return list(self._model_loaders.keys())
    
    def get_feat_map(self, name: str) -> Optional[str]:
        """Get the featurizer mapping for a model."""
        # Trigger model load if not already loaded
        if name not in self._loaded_models and name in self._model_loaders:
            try:
                self[name]  # This loads the model and updates feat map
            except ImportError:
                pass
        return self._model_feat_map.get(name)


# Global lazy registry instance
model_address_map = LazyModelRegistry()

# For backward compatibility
MODEL_FEAT_MAP = model_address_map._model_feat_map
model_names = model_address_map.keys()

