"""
Primitives client module for interacting with DeepChem Server API endpoints.

Contains the Primitives class for all computation tasks.
"""

from typing import Any, Dict, List, Optional, Union

from .base import BaseClient
from .settings import Settings

class Primitives(BaseClient):
    """
    Client for interacting with DeepChem Server primitive endpoints.

    This class provides methods to submit jobs to various primitive endpoints
    including featurization, training, evaluation, and inference.
    """

    def __init__(self, settings: Optional[Settings] = None, base_url: Optional[str] = None):
        """
        Initialize Primitives client.

        Args:
            settings: Settings instance for configuration
            base_url: Base URL for the API (overrides settings if provided)
        """
        super().__init__(settings, base_url)

    def featurize(
        self,
        dataset_address: str,
        featurizer: str,
        output: str,
        dataset_column: str,
        feat_kwargs: Optional[Dict[str, Any]] = None,
        label_column: Optional[str] = None,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Submit a featurization job.

        Args:
            dataset_address: Datastore address of dataset to featurize
            featurizer: Featurizer to use
            output: Name of the featurized dataset
            dataset_column: Column containing the input for featurizer
            feat_kwargs: Keyword arguments to pass to featurizer on initialization
            label_column: The target column in case this dataset is going to be used for training
            profile_name: Profile name (uses settings if not provided)
            project_name: Project name (uses settings if not provided)

        Returns:
            Response containing the featurized file address

        Raises:
            ValueError: If required settings are missing
            requests.exceptions.RequestException: If API request fails
        """
        # Get profile and project names (validates configuration)
        profile, project = self._get_profile_and_project(profile_name, project_name)

        if feat_kwargs is None:
            feat_kwargs = {}

        # Prepare the request data
        data = {
            "profile_name": profile,
            "project_name": project,
            "dataset_address": dataset_address,
            "featurizer": featurizer,
            "output": output,
            "dataset_column": dataset_column,
            "feat_kwargs": {
                "feat_kwargs": feat_kwargs
            },
        }

        if label_column is not None:
            data["label_column"] = label_column

        response = self._post("/primitive/featurize", json=data)
        return self._validate_response(response)

    def train(
        self,
        dataset_address: str,
        model_type: str,
        model_name: str,
        init_kwargs: Optional[Dict[str, Any]] = None,
        train_kwargs: Optional[Dict[str, Any]] = None,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Submit a training job.

        Args:
            dataset_address: Datastore address of dataset to train on
            model_type: Type of model to train
            model_name: Name of the trained model
            init_kwargs: Keyword arguments to pass to model on initialization
            train_kwargs: Keyword arguments to pass to model on training
            profile_name: Profile name (uses settings if not provided)
            project_name: Project name (uses settings if not provided)

        Returns:
            Response containing the trained model address

        Raises:
            ValueError: If required settings are missing
            requests.exceptions.RequestException: If API request fails
        """
        # Get profile and project names (validates configuration)
        profile, project = self._get_profile_and_project(profile_name, project_name)

        if init_kwargs is None:
            init_kwargs = {}
        if train_kwargs is None:
            train_kwargs = {}

        # Prepare the request data
        data = {
            "profile_name": profile,
            "project_name": project,
            "dataset_address": dataset_address,
            "model_type": model_type,
            "model_name": model_name,
            "init_kwargs": init_kwargs,
            "train_kwargs": train_kwargs,
        }

        response = self._post("/primitive/train", json=data)
        return self._validate_response(response)

    def evaluate(
        self,
        dataset_addresses: List[str],
        model_address: str,
        metrics: List[str],
        output_key: str,
        is_metric_plots: bool = False,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Submit an evaluation job.

        Args:
            dataset_addresses: List of dataset addresses to evaluate the model on
            model_address: Datastore address of the trained model
            metrics: List of metrics to evaluate the model with
            output_key: Name of the evaluation output
            is_metric_plots: Whether plot based metric is used or not
            profile_name: Profile name (uses settings if not provided)
            project_name: Project name (uses settings if not provided)

        Returns:
            Response containing the evaluation result address

        Raises:
            ValueError: If required settings are missing
            requests.exceptions.RequestException: If API request fails
        """
        # Get profile and project names (validates configuration)
        profile, project = self._get_profile_and_project(profile_name, project_name)

        # Prepare the request data
        data = {
            "profile_name": profile,
            "project_name": project,
            "dataset_addresses": dataset_addresses,
            "model_address": model_address,
            "metrics": metrics,
            "output_key": output_key,
            "is_metric_plots": is_metric_plots,
        }

        response = self._post("/primitive/evaluate", json=data)
        return self._validate_response(response)

    def infer(
        self,
        model_address: str,
        data_address: str,
        output: str,
        dataset_column: Optional[str] = None,
        shard_size: Optional[int] = 8192,
        threshold: Optional[Union[int, float]] = None,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Submit an inference job.

        Args:
            model_address: Datastore address of the trained model
            data_address: Datastore address of the dataset to run inference on
            output: Name of the output inference results
            dataset_column: Column in the dataset to perform inference on (required for raw CSV data)
            shard_size: Shard size for the inference operation (default: 8192)
            threshold: Threshold for binarizing predictions (for classification tasks)
            profile_name: Profile name (uses settings if not provided)
            project_name: Project name (uses settings if not provided)

        Returns:
            Response containing the inference results address

        Raises:
            ValueError: If required settings are missing
            requests.exceptions.RequestException: If API request fails
        """
        # Get profile and project names (validates configuration)
        profile, project = self._get_profile_and_project(profile_name, project_name)

        # Prepare the request data
        data = {
            "profile_name": profile,
            "project_name": project,
            "model_address": model_address,
            "data_address": data_address,
            "output": output,
            "dataset_column": dataset_column,
            "shard_size": shard_size,
            "threshold": threshold,
        }

        response = self._post("/primitive/infer", json=data)
        return self._validate_response(response)
