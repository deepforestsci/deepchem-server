"""
Tests for the Primitives module.
"""

from unittest.mock import Mock, patch

import pytest
import responses

from pyds import Settings
from pyds.primitives import Primitives


class TestPrimitivesCommonUtils:
    """Test cases for Primitives class common utilities and initialization."""

    def test_init_with_settings(self, settings_instance):
        """Test Primitives initialization with settings."""
        primitives_client = Primitives(settings=settings_instance)

        assert primitives_client.settings == settings_instance
        assert primitives_client.base_url == "http://localhost:8000"

    def test_init_with_base_url_override(self, settings_instance):
        """Test Primitives initialization with base_url override."""
        override_url = "http://override:9000"
        primitives_client = Primitives(settings=settings_instance,
                                       base_url=override_url)

        assert primitives_client.base_url == override_url

    def test_init_without_settings(self):
        """Test Primitives initialization without settings."""
        with patch("pyds.base.client.Settings") as mock_settings:
            mock_instance = Mock()
            mock_instance.get_base_url.return_value = "http://default:8000"
            mock_settings.return_value = mock_instance

            primitives_client = Primitives()

            assert primitives_client.settings == mock_instance

    @responses.activate
    def test_api_error_handling(self, primitives_client_instance):
        """Test API error handling across all methods."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/featurize",
            json={"detail": "Featurizer not found"},
            status=400,
        )

        with pytest.raises(Exception) as exc_info:
            primitives_client_instance.featurize(
                dataset_address="deepchem://test-profile/test-project/data.csv",
                featurizer="invalid_featurizer",
                output="featurized_data",
                dataset_column="smiles",
            )

        assert "Featurizer not found" in str(exc_info.value)

    @responses.activate
    def test_server_error_handling(self, primitives_client_instance):
        """Test server error handling."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/train",
            json={"detail": "Internal server error"},
            status=500,
        )

        with pytest.raises(Exception) as exc_info:
            primitives_client_instance.train(
                dataset_address="deepchem://test-profile/test-project/data.csv",
                model_type="multitask_classifier",
                model_name="test_model",
            )

        assert "Internal server error" in str(exc_info.value)

    def test_all_methods_require_configuration(self, temp_settings_file):
        """Test that all methods require proper configuration."""
        settings = Settings(settings_file=temp_settings_file)
        primitives_client = Primitives(settings=settings)

        methods_to_test = [
            (
                "featurize",
                {
                    "dataset_address": "test",
                    "featurizer": "ecfp",
                    "output": "test",
                    "dataset_column": "smiles",
                },
            ),
            (
                "train",
                {
                    "dataset_address": "test",
                    "model_type": "classifier",
                    "model_name": "test"
                },
            ),
            (
                "evaluate",
                {
                    "dataset_addresses": ["test"],
                    "model_address": "test",
                    "metrics": ["accuracy"],
                    "output_key": "test",
                },
            ),
            ("infer", {
                "model_address": "test",
                "data_address": "test",
                "output": "test"
            }),
        ]

        for method_name, kwargs in methods_to_test:
            with pytest.raises(ValueError) as exc_info:
                method = getattr(primitives_client, method_name)
                method(**kwargs)

            assert "Missing required settings" in str(exc_info.value)


class TestFeaturize:
    """Test cases for the featurize primitive operation."""

    @responses.activate
    def test_featurize_success(self, primitives_client_instance):
        """Test successful featurization."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/featurize",
            json={
                "status":
                    "success",
                "featurized_file_address":
                    "deepchem://test-profile/test-project/featurized_data",
            },
            status=200,
        )

        result = primitives_client_instance.featurize(
            dataset_address="deepchem://test-profile/test-project/data.csv",
            featurizer="ecfp",
            output="featurized_data",
            dataset_column="smiles",
            label_column="logp",
            feat_kwargs={"size": 1024},
        )

        assert result["status"] == "success"
        assert "featurized_file_address" in result

        # Verify request payload
        request_body = responses.calls[0].request.body
        assert b'"featurizer": "ecfp"' in request_body
        assert b'"dataset_column": "smiles"' in request_body

    @responses.activate
    def test_featurize_with_custom_profile_project(self,
                                                   primitives_client_instance):
        """Test featurization with custom profile and project."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/featurize",
            json={"status": "success"},
            status=200,
        )

        primitives_client_instance.featurize(
            dataset_address="deepchem://custom-profile/custom-project/data.csv",
            featurizer="mordred",
            output="featurized_data",
            dataset_column="smiles",
            profile_name="custom-profile",
            project_name="custom-project",
        )

        request_body = responses.calls[0].request.body.decode("utf-8")
        assert '"profile_name": "custom-profile"' in request_body
        assert '"project_name": "custom-project"' in request_body

    @responses.activate
    def test_featurize_without_label_column(self, primitives_client_instance):
        """Test featurization without label column."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/featurize",
            json={"status": "success"},
            status=200,
        )

        primitives_client_instance.featurize(
            dataset_address="deepchem://test-profile/test-project/data.csv",
            featurizer="ecfp",
            output="featurized_data",
            dataset_column="smiles",
        )

        request_body = responses.calls[0].request.body.decode("utf-8")
        assert '"label_column"' not in request_body

    @responses.activate
    def test_featurize_with_empty_feat_kwargs(self, primitives_client_instance):
        """Test featurization with empty feat_kwargs."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/featurize",
            json={"status": "success"},
            status=200,
        )

        primitives_client_instance.featurize(
            dataset_address="deepchem://test-profile/test-project/data.csv",
            featurizer="ecfp",
            output="featurized_data",
            dataset_column="smiles",
        )

        request_body = responses.calls[0].request.body.decode("utf-8")
        assert '"feat_kwargs": {"feat_kwargs": {}}' in request_body

    def test_featurize_requires_configuration(self, temp_settings_file):
        """Test that featurize requires proper configuration."""
        settings = Settings(settings_file=temp_settings_file)
        primitives_client = Primitives(settings=settings)

        with pytest.raises(ValueError) as exc_info:
            primitives_client.featurize(
                dataset_address="deepchem://test/test/data.csv",
                featurizer="ecfp",
                output="featurized_data",
                dataset_column="smiles",
            )

        assert "Missing required settings" in str(exc_info.value)


class TestTrain:
    """Test cases for the train primitive operation."""

    @responses.activate
    def test_train_success(self, primitives_client_instance):
        """Test successful training."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/train",
            json={
                "status":
                    "success",
                "model_address":
                    "deepchem://test-profile/test-project/trained_model",
            },
            status=200,
        )

        dataset_address = "deepchem://test-profile/test-project/featurized_data"
        result = primitives_client_instance.train(
            dataset_address=dataset_address,
            model_type="multitask_classifier",
            model_name="trained_model",
            init_kwargs={
                "n_tasks": 1,
                "n_features": 1024
            },
            train_kwargs={"nb_epoch": 10},
        )

        assert result["status"] == "success"
        assert "model_address" in result

        # Verify request payload
        request_body = responses.calls[0].request.body.decode("utf-8")
        assert '"model_type": "multitask_classifier"' in request_body
        assert '"model_name": "trained_model"' in request_body

    @responses.activate
    def test_train_with_empty_kwargs(self, primitives_client_instance):
        """Test training with empty init_kwargs and train_kwargs."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/train",
            json={"status": "success"},
            status=200,
        )

        dataset_address = "deepchem://test-profile/test-project/featurized_data"
        primitives_client_instance.train(
            dataset_address=dataset_address,
            model_type="multitask_classifier",
            model_name="trained_model",
        )

        request_body = responses.calls[0].request.body.decode("utf-8")
        assert '"init_kwargs": {}' in request_body
        assert '"train_kwargs": {}' in request_body


class TestEvaluate:
    """Test cases for the evaluate primitive operation."""

    @responses.activate
    def test_evaluate_success(self, primitives_client_instance):
        """Test successful evaluation."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/evaluate",
            json={
                "status":
                    "success",
                "evaluation_address":
                    "deepchem://test-profile/test-project/evaluation_results",
            },
            status=200,
        )

        result = primitives_client_instance.evaluate(
            dataset_addresses=[
                "deepchem://test-profile/test-project/test_data"
            ],
            model_address="deepchem://test-profile/test-project/trained_model",
            metrics=["roc_auc_score", "accuracy_score"],
            output_key="evaluation_results",
            is_metric_plots=True,
        )

        assert result["status"] == "success"
        assert "evaluation_address" in result

        # Verify request payload
        request_body = responses.calls[0].request.body.decode("utf-8")
        assert '"metrics": ["roc_auc_score", "accuracy_score"]' in request_body
        assert '"is_metric_plots": true' in request_body

    @responses.activate
    def test_evaluate_with_multiple_datasets(self, primitives_client_instance):
        """Test evaluation with multiple datasets."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/evaluate",
            json={"status": "success"},
            status=200,
        )

        dataset_addresses = [
            "deepchem://test-profile/test-project/train_data",
            "deepchem://test-profile/test-project/test_data",
            "deepchem://test-profile/test-project/valid_data",
        ]

        primitives_client_instance.evaluate(
            dataset_addresses=dataset_addresses,
            model_address="deepchem://test-profile/test-project/trained_model",
            metrics=["accuracy_score"],
            output_key="evaluation_results",
        )

        request_body = responses.calls[0].request.body.decode("utf-8")
        assert (
            '"dataset_addresses": ["deepchem://test-profile/test-project/train_data"'
            in request_body)

    @responses.activate
    def test_evaluate_default_is_metric_plots(self, primitives_client_instance):
        """Test evaluation with default is_metric_plots value."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/evaluate",
            json={"status": "success"},
            status=200,
        )

        primitives_client_instance.evaluate(
            dataset_addresses=[
                "deepchem://test-profile/test-project/test_data"
            ],
            model_address="deepchem://test-profile/test-project/trained_model",
            metrics=["accuracy_score"],
            output_key="evaluation_results",
        )

        request_body = responses.calls[0].request.body.decode("utf-8")
        assert '"is_metric_plots": false' in request_body


class TestInfer:
    """Test cases for the infer primitive operation."""

    @responses.activate
    def test_infer_success(self, primitives_client_instance):
        """Test successful inference."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/infer",
            json={
                "status":
                    "success",
                "inference_address":
                    "deepchem://test-profile/test-project/inference_results",
            },
            status=200,
        )

        result = primitives_client_instance.infer(
            model_address="deepchem://test-profile/test-project/trained_model",
            data_address="deepchem://test-profile/test-project/new_data",
            output="inference_results",
            dataset_column="smiles",
            shard_size=4096,
            threshold=0.5,
        )

        assert result["status"] == "success"
        assert "inference_address" in result

        # Verify request payload
        request_body = responses.calls[0].request.body.decode("utf-8")
        assert '"shard_size": 4096' in request_body
        assert '"threshold": 0.5' in request_body

    @responses.activate
    def test_infer_with_default_values(self, primitives_client_instance):
        """Test inference with default parameter values."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/infer",
            json={"status": "success"},
            status=200,
        )

        primitives_client_instance.infer(
            model_address="deepchem://test-profile/test-project/trained_model",
            data_address="deepchem://test-profile/test-project/new_data",
            output="inference_results",
        )

        request_body = responses.calls[0].request.body.decode("utf-8")
        assert '"shard_size": 8192' in request_body  # Default value
        assert '"dataset_column": null' in request_body
        assert '"threshold": null' in request_body

    @responses.activate
    def test_infer_without_optional_params(self, primitives_client_instance):
        """Test inference without optional parameters."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/infer",
            json={"status": "success"},
            status=200,
        )

        primitives_client_instance.infer(
            model_address="deepchem://test-profile/test-project/trained_model",
            data_address="deepchem://test-profile/test-project/new_data",
            output="inference_results",
            dataset_column="smiles",
        )

        request_body = responses.calls[0].request.body.decode("utf-8")
        assert '"dataset_column": "smiles"' in request_body
