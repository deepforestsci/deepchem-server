"""
Integration tests for upload and featurize functionality.

These tests require a running DeepChem server and are intended for integration testing.
They have been refactored to use pytest and the pyds client library.
"""

import os

import pytest
import responses

from pyds import Data, Primitives, Settings


class TestIntegrationUploadFeaturize:
    """Integration tests for upload and featurize workflow."""

    @pytest.fixture
    def integration_settings(self, temp_settings_file):
        """Create settings for integration testing."""
        return Settings(
            settings_file=temp_settings_file,
            profile="test-profile",
            project="test-project",
            base_url="http://0.0.0.0:8000",
        )

    @pytest.fixture
    def data_client(self, integration_settings):
        """Create Data client for integration testing."""
        return Data(settings=integration_settings)

    @pytest.fixture
    def primitives_client(self, integration_settings):
        """Create Primitives client for integration testing."""
        return Primitives(settings=integration_settings)

    @pytest.fixture
    def zinc10_csv_file(self, test_assets_dir):
        """Path to the zinc10.csv test file."""
        return test_assets_dir / "zinc10.csv"

    @pytest.mark.integration
    def test_upload_csv_integration(self, data_client, zinc10_csv_file):
        """
        Integration test for CSV file upload.

        Note: This test requires a running DeepChem server.
        Mark with @pytest.mark.integration to run separately.
        """
        if not zinc10_csv_file.exists():
            pytest.skip("zinc10.csv test file not found")

        result = data_client.upload_data(
            file_path=zinc10_csv_file,
            filename="zinc10_sample.csv",
            description="Sample test csv file",
        )

        expected_address = "deepchem://test-profile/test-project/zinc10_sample.csv"
        assert result.get("dataset_address") == expected_address

    @pytest.mark.integration
    def test_featurize_integration(self, data_client, primitives_client,
                                   zinc10_csv_file):
        """
        Integration test for featurization workflow.

        Note: This test requires a running DeepChem server.
        Mark with @pytest.mark.integration to run separately.
        """
        if not zinc10_csv_file.exists():
            pytest.skip("zinc10.csv test file not found")

        # First upload the file
        upload_result = data_client.upload_data(
            file_path=zinc10_csv_file,
            filename="zinc10_sample.csv",
            description="Sample test csv file",
        )

        dataset_address = upload_result["dataset_address"]

        # Then featurize it
        featurize_result = primitives_client.featurize(
            dataset_address=dataset_address,
            featurizer="ecfp",
            output="test_featurized",
            dataset_column="smiles",
            label_column="logp",
            feat_kwargs={"size": 1024},
        )

        expected_address = "deepchem://test-profile/test-project/test_featurized"
        assert featurize_result.get(
            "featurized_file_address") == expected_address

    @responses.activate
    def test_upload_csv_mocked(self, data_client, zinc10_csv_file):
        """
        Mocked test for CSV file upload functionality.

        This test doesn't require a running server and uses mocked responses.
        """
        # Mock the upload response
        responses.add(
            responses.POST,
            "http://0.0.0.0:8000/data/uploaddata",
            json={
                "status":
                    "success",
                "dataset_address":
                    "deepchem://test-profile/test-project/zinc10_sample.csv",
            },
            status=200,
        )

        # Create a temporary CSV file for testing if zinc10.csv doesn't exist
        if not zinc10_csv_file.exists():
            import tempfile

            content = "smiles,logp\nCCO,0.2\nCCC,1.0\n"
            with tempfile.NamedTemporaryFile(mode="w",
                                             suffix=".csv",
                                             delete=False) as f:
                f.write(content)
                test_file = f.name
        else:
            test_file = str(zinc10_csv_file)

        try:
            result = data_client.upload_data(
                file_path=test_file,
                filename="zinc10_sample.csv",
                description="Sample test csv file",
            )

            assert result["status"] == "success"
            assert (result["dataset_address"] ==
                    "deepchem://test-profile/test-project/zinc10_sample.csv")
        finally:
            # Cleanup temporary file if created
            if not zinc10_csv_file.exists() and os.path.exists(test_file):
                os.unlink(test_file)

    @responses.activate
    def test_featurize_mocked(self, primitives_client):
        """
        Mocked test for featurization functionality.

        This test doesn't require a running server and uses mocked responses.
        """
        # Mock the featurize response
        responses.add(
            responses.POST,
            "http://0.0.0.0:8000/primitive/featurize",
            json={
                "status":
                    "success",
                "featurized_file_address":
                    "deepchem://test-profile/test-project/test_featurized",
            },
            status=200,
        )

        result = primitives_client.featurize(
            dataset_address=
            "deepchem://test-profile/test-project/zinc10_sample.csv",
            featurizer="ecfp",
            output="test_featurized",
            dataset_column="smiles",
            label_column="logp",
            feat_kwargs={"size": 1024},
        )

        assert result["status"] == "success"
        assert (result["featurized_file_address"] ==
                "deepchem://test-profile/test-project/test_featurized")

    @responses.activate
    def test_upload_error_handling(self, data_client, sample_csv_file):
        """Test error handling in upload functionality."""
        # Mock an error response
        responses.add(
            responses.POST,
            "http://0.0.0.0:8000/data/uploaddata",
            json={"detail": "File format not supported"},
            status=400,
        )

        with pytest.raises(Exception) as exc_info:
            data_client.upload_data(file_path=sample_csv_file,
                                    filename="test.csv")

        assert "File format not supported" in str(exc_info.value)

    @responses.activate
    def test_featurize_error_handling(self, primitives_client):
        """Test error handling in featurize functionality."""
        # Mock an error response
        responses.add(
            responses.POST,
            "http://0.0.0.0:8000/primitive/featurize",
            json={"detail": "Invalid featurizer"},
            status=400,
        )

        with pytest.raises(Exception) as exc_info:
            primitives_client.featurize(
                dataset_address="deepchem://test-profile/test-project/data.csv",
                featurizer="invalid_featurizer",
                output="test_output",
                dataset_column="smiles",
            )

        assert "Invalid featurizer" in str(exc_info.value)
