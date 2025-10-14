"""
Unit tests for Train primitive using live server.
"""

import time
from typing import Any
import uuid

import pytest

from pyds.data import Data
from pyds.primitives import Train
from pyds.settings import Settings


class TestTrain:
    """Unit tests for Train primitive."""

    def test_init(self, test_settings: Settings) -> None:
        """Test Train initialization."""
        client = Train(settings=test_settings)

        assert client.settings == test_settings
        assert client.base_url == "http://localhost:8000"

    def test_train_random_forest_classifier(
        self,
        live_train_client: Train,
        live_data_client: Data,
        live_featurize_client: Any,
        large_classification_csv: str,
    ) -> None:
        """Test training a random forest classifier with real data."""
        # Generate unique identifiers to avoid naming conflicts
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        test_file = large_classification_csv

        upload_result = live_data_client.upload_data(
            file_path=test_file,
            filename=f"train_rf_cls_{test_id}_{timestamp}.csv",
            description="Test data for random forest classifier training",
        )

        # Featurize data using ECFP with unique output name
        featurize_result = live_featurize_client.run(
            dataset_address=upload_result["dataset_address"],
            featurizer="ecfp",
            output=f"train_rf_feat_{test_id}_{timestamp}",
            dataset_column="smiles",
            label_column="label",
            feat_kwargs={
                "radius": 2,
                "size": 1024
            },
        )

        result = live_train_client.run(
            dataset_address=featurize_result["featurized_file_address"],
            model_type="random_forest_classifier",
            model_name=f"rf_cls_model_{test_id}_{timestamp}",
            init_kwargs={
                "n_estimators": 50,
                "random_state": 42
            },
            train_kwargs={},
        )

        assert "trained_model_address" in result
        assert isinstance(result["trained_model_address"], str)
        assert len(result["trained_model_address"]) > 0

    def test_train_random_forest_regressor(
        self,
        live_train_client: Train,
        live_data_client: Data,
        live_featurize_client: Any,
        complex_regression_csv: str,
    ) -> None:
        """Test training a random forest regressor with real data."""
        # Generate unique identifiers to avoid naming conflicts
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        test_file = complex_regression_csv

        upload_result = live_data_client.upload_data(
            file_path=test_file,
            filename=f"train_rf_reg_{test_id}_{timestamp}.csv",
            description="Test data for random forest regressor training",
        )

        featurize_result = live_featurize_client.run(
            dataset_address=upload_result["dataset_address"],
            featurizer="ecfp",
            output=f"train_rf_reg_feat_{test_id}_{timestamp}",
            dataset_column="smiles",
            label_column="target",
            feat_kwargs={
                "radius": 2,
                "size": 512
            },
        )

        result = live_train_client.run(
            dataset_address=featurize_result["featurized_file_address"],
            model_type="random_forest_regressor",
            model_name=f"rf_reg_model_{test_id}_{timestamp}",
            init_kwargs={
                "n_estimators": 30,
                "random_state": 42
            },
            train_kwargs={},
        )

        # Verify the response structure
        assert "trained_model_address" in result
        assert isinstance(result["trained_model_address"], str)
        assert len(result["trained_model_address"]) > 0

    def test_train_linear_regression(
        self,
        live_train_client: Train,
        live_data_client: Data,
        live_featurize_client: Any,
        simple_regression_csv: str,
    ) -> None:
        """Test training a linear regression model with real data."""
        # Generate unique identifiers to avoid naming conflicts
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        test_file = simple_regression_csv

        upload_result = live_data_client.upload_data(
            file_path=test_file,
            filename=f"train_linear_reg_{test_id}_{timestamp}.csv",
            description="Test data for linear regression training",
        )

        featurize_result = live_featurize_client.run(
            dataset_address=upload_result["dataset_address"],
            featurizer="ecfp",
            output=f"train_linear_feat_{test_id}_{timestamp}",
            dataset_column="smiles",
            label_column="property",
        )

        result = live_train_client.run(
            dataset_address=featurize_result["featurized_file_address"],
            model_type="linear_regression",
            model_name=f"linear_model_{test_id}_{timestamp}",
            init_kwargs={"fit_intercept": True},
            train_kwargs={},
        )

        # Verify the response structure
        assert "trained_model_address" in result
        assert isinstance(result["trained_model_address"], str)
        assert len(result["trained_model_address"]) > 0

    def test_train_with_minimal_parameters(
        self,
        live_train_client: Train,
        live_data_client: Data,
        live_featurize_client: Any,
        minimal_classification_csv: str,
    ) -> None:
        """Test training with minimal parameters (using defaults)."""
        # Generate unique identifiers to avoid naming conflicts
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        test_file = minimal_classification_csv

        upload_result = live_data_client.upload_data(
            file_path=test_file,
            filename=f"train_minimal_{test_id}_{timestamp}.csv",
        )

        featurize_result = live_featurize_client.run(
            dataset_address=upload_result["dataset_address"],
            featurizer="ecfp",
            output=f"train_minimal_feat_{test_id}_{timestamp}",
            dataset_column="smiles",
            label_column="label",
        )

        # Train with minimal parameters (no init_kwargs, no train_kwargs)
        result = live_train_client.run(
            dataset_address=featurize_result["featurized_file_address"],
            model_type="random_forest_classifier",
            model_name=f"minimal_model_{test_id}_{timestamp}",
        )

        assert "trained_model_address" in result
        assert isinstance(result["trained_model_address"], str)

    def test_train_invalid_model_type(self, live_train_client: Train) -> None:
        """Test training with invalid model type."""
        with pytest.raises(Exception, match="Invalid model type"):
            live_train_client.run(
                dataset_address="test/dataset",
                model_type="NonExistentModel",
                model_name="test_model",
            )

    def test_train_nonexistent_dataset(self, live_train_client: Train) -> None:
        """Test training with non-existent dataset."""
        with pytest.raises(Exception):
            live_train_client.run(
                dataset_address="nonexistent/dataset/address",
                model_type="random_forest_classifier",
                model_name="test_model",
            )

    def test_run_missing_settings(self, test_settings_not_configured: Settings) -> None:
        """Test train run with missing settings."""
        client = Train(settings=test_settings_not_configured)

        with pytest.raises(ValueError, match="Missing required settings"):
            client.run(
                dataset_address="test/dataset",
                model_type="random_forest_classifier",
                model_name="trained_model",
            )
