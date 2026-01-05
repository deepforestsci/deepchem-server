"""
Unit tests for Featurize primitive using live server.
"""

import time
import uuid

import pytest

from pyds.data import Data
from pyds.models import DeepchemData, Job
from pyds.primitives import Featurize
from pyds.settings import Settings


class TestFeaturize:
    """Unit tests for Featurize primitive."""

    def test_init(self, test_settings: Settings) -> None:
        """Test Featurize initialization."""
        client = Featurize(settings=test_settings)

        assert client.settings == test_settings
        assert client.base_url == "http://localhost:8000"

    def test_run_success(
        self,
        live_featurize_client: Featurize,
        live_data_client: Data,
        small_classification_csv: str,
    ) -> None:
        """Test successful featurize run on live server."""
        # Generate unique identifiers to avoid naming conflicts
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        test_file = small_classification_csv

        upload_result = live_data_client.upload_data(
            file_path=test_file,
            filename=f"test_featurize_{test_id}_{timestamp}.csv",
            description="Test data for featurization",
        )
        # Use .address instead of ["dataset_address"]
        dataset_address = upload_result.address

        # Test featurization - returns Job object
        job = live_featurize_client.run(
            dataset_address=dataset_address,
            featurizer="ecfp",
            output=f"test_featurized_output_{test_id}_{timestamp}",
            dataset_column="smiles",
            feat_kwargs={
                "radius": 2,
                "size": 1024
            },
            label_column="label",
        )

        # Verify Job object returned
        assert isinstance(job, Job)
        assert job.id is not None and job.id != ""

        # wait for job to complete
        job.wait()
        featurized_data = job.result_as_data()
        assert isinstance(featurized_data, DeepchemData)
        assert featurized_data.address is not None and featurized_data.address != ""
        assert f"test_featurized_output_{test_id}_{timestamp}" in featurized_data.address

    def test_run_with_defaults(
        self,
        live_featurize_client: Featurize,
        live_data_client: Data,
        small_classification_csv: str,
    ) -> None:
        """Test featurize run with default parameters on live server."""
        # Generate unique identifiers to avoid naming conflicts
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        test_file = small_classification_csv

        upload_result = live_data_client.upload_data(
            file_path=test_file,
            filename=f"test_featurize_defaults_{test_id}_{timestamp}.csv",
            description="Test data for featurization with defaults",
        )
        # Use .address instead of ["dataset_address"]
        dataset_address = upload_result.address

        # Test featurization with minimal parameters - returns Job object
        job = live_featurize_client.run(
            dataset_address=dataset_address,
            featurizer="ecfp",
            output=f"test_featurized_defaults_{test_id}_{timestamp}",
            dataset_column="smiles",
        )

        # Verify Job object returned
        assert isinstance(job, Job)
        assert job.id is not None and job.id != ""

        # wait for job to complete
        job.wait()
        featurized_data = job.result_as_data()
        assert isinstance(featurized_data, DeepchemData)
        assert featurized_data.address is not None and featurized_data.address != ""
        assert f"test_featurized_defaults_{test_id}_{timestamp}" in featurized_data.address

    def test_run_with_profile_project_override(
        self,
        live_featurize_client: Featurize,
        live_data_client: Data,
        small_classification_csv: str,
    ) -> None:
        """Test featurize run with profile and project override on live server."""
        # Generate unique identifiers to avoid naming conflicts
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        test_file = small_classification_csv

        upload_result = live_data_client.upload_data(
            file_path=test_file,
            filename=f"test_featurize_override_{test_id}_{timestamp}.csv",
            description="Test data for featurization with override",
        )
        # Use .address instead of ["dataset_address"]
        dataset_address = upload_result.address

        # Test featurization with profile/project override - returns Job object
        job = live_featurize_client.run(
            dataset_address=dataset_address,
            featurizer="ecfp",
            output=f"test_featurized_override_{test_id}_{timestamp}",
            dataset_column="smiles",
            profile_name="test_profile",
            project_name="test_project",
        )

        # Verify Job object returned
        assert isinstance(job, Job)
        assert job.id is not None and job.id != ""

        # wait for job to complete
        job.wait()
        featurized_data = job.result_as_data()
        assert isinstance(featurized_data, DeepchemData)
        assert featurized_data.address is not None and featurized_data.address != ""
        assert f"test_featurized_override_{test_id}_{timestamp}" in featurized_data.address

    def test_run_missing_settings(self, test_settings_not_configured: Settings) -> None:
        """Test featurize run with missing settings."""
        client = Featurize(settings=test_settings_not_configured)

        with pytest.raises(ValueError, match="Missing required settings"):
            client.run(
                dataset_address="test/dataset",
                featurizer="CircularFingerprint",
                output="featurized_output",
                dataset_column="smiles",
            )

    def test_run_api_error(self, live_featurize_client: Featurize) -> None:
        """Test featurize run with API error on live server."""
        # Note: Server now queues job even with invalid featurizer
        # The job will fail during execution, so we just verify a Job is returned
        job = live_featurize_client.run(
            dataset_address="test/dataset",
            featurizer="InvalidFeaturizer",
            output="featurized_output",
            dataset_column="smiles",
        )
        # Job is created but will fail during worker execution
        assert isinstance(job, Job)
        assert job.id is not None and job.id != ""

        # wait for job to fail - should raise exception with error message
        with pytest.raises(Exception, match="Featurizer not recognized"):
            job.wait()
