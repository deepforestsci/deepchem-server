"""
Primitives client module for interacting with DeepChem server API endpoints.

Contains the Primitives class for all computation tasks.
"""

import json
from typing import Optional, Dict, Any
from .base import BaseClient
from .settings import Settings


class Primitives(BaseClient):
    """
    Client for interacting with DeepChem server primitive endpoints.

    This class provides methods to submit jobs to various primitive endpoints
    such as featurization.
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
            "feat_kwargs": {"feat_kwargs": feat_kwargs},
        }

        if label_column is not None:
            data["label_column"] = label_column

        return self._make_request("POST", "/primitive/featurize", json=data)
