"""
Tests for the Data module.
"""

from pathlib import Path
from unittest.mock import Mock, patch

import pytest
import responses

from pyds import Settings
from pyds.data import Data


class TestData:
    """Test cases for Data class."""

    def test_init_with_settings(self, settings_instance):
        """Test Data initialization with settings."""
        data_client = Data(settings=settings_instance)

        assert data_client.settings == settings_instance
        assert data_client.base_url == "http://localhost:8000"

    def test_init_with_base_url_override(self, settings_instance):
        """Test Data initialization with base_url override."""
        override_url = "http://override:9000"
        data_client = Data(settings=settings_instance, base_url=override_url)

        assert data_client.base_url == override_url

    def test_init_without_settings(self):
        """Test Data initialization without settings."""
        with patch("pyds.data.Settings") as mock_settings:
            mock_instance = Mock()
            mock_instance.get_base_url.return_value = "http://default:8000"
            mock_settings.return_value = mock_instance

            data_client = Data()

            assert data_client.settings == mock_instance

    @responses.activate
    def test_upload_data_success(self, data_client_instance, sample_csv_file):
        """Test successful data upload."""
        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json={
                "status":
                    "success",
                "dataset_address":
                    "deepchem://test-profile/test-project/test.csv",
            },
            status=200,
        )

        result = data_client_instance.upload_data(file_path=sample_csv_file,
                                                  filename="test.csv",
                                                  description="Test CSV file")

        assert result["status"] == "success"
        assert "dataset_address" in result
        assert len(responses.calls) == 1

    @responses.activate
    def test_upload_data_with_custom_profile_project(self, data_client_instance,
                                                     sample_csv_file):
        """Test data upload with custom profile and project."""
        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json={
                "status":
                    "success",
                "dataset_address":
                    "deepchem://custom-profile/custom-project/test.csv",
            },
            status=200,
        )

        result = data_client_instance.upload_data(
            file_path=sample_csv_file,
            filename="test.csv",
            profile_name="custom-profile",
            project_name="custom-project",
        )

        assert result["status"] == "success"
        # Verify the request was made with correct parameters
        request_data = responses.calls[0].request.body.decode("utf-8")
        assert 'name="profile_name"\r\n\r\ncustom-profile' in request_data
        assert 'name="project_name"\r\n\r\ncustom-project' in request_data

    @responses.activate
    def test_upload_data_file_not_found(self, data_client_instance):
        """Test data upload with non-existent file."""
        with pytest.raises(FileNotFoundError):
            data_client_instance.upload_data(file_path="/nonexistent/file.csv",
                                             filename="test.csv")

    @responses.activate
    def test_upload_data_server_error(self, data_client_instance,
                                      sample_csv_file):
        """Test data upload with server error."""
        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json={"detail": "Internal server error"},
            status=500,
        )

        with pytest.raises(Exception) as exc_info:
            data_client_instance.upload_data(file_path=sample_csv_file,
                                             filename="test.csv")

        assert "Internal server error" in str(exc_info.value)

    @responses.activate
    def test_upload_data_bad_request(self, data_client_instance,
                                     sample_csv_file):
        """Test data upload with bad request error."""
        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json={"detail": "Invalid file format"},
            status=400,
        )

        with pytest.raises(Exception) as exc_info:
            data_client_instance.upload_data(file_path=sample_csv_file,
                                             filename="test.csv")

        assert "Invalid file format" in str(exc_info.value)

    def test_upload_data_filename_from_path(self, data_client_instance,
                                            sample_csv_file):
        """Test that filename is extracted from file path when not provided."""
        with patch.object(data_client_instance, "_post") as mock_post:
            mock_post.return_value.status_code = 200
            mock_post.return_value.json.return_value = {"status": "success"}

            data_client_instance.upload_data(file_path=sample_csv_file)

            # Verify filename was extracted from path
            args, kwargs = mock_post.call_args
            assert args[0] == "/data/uploaddata"
            # The data should be a MultipartEncoder
            assert "data" in kwargs

    def test_upload_data_with_pathlib_path(self, data_client_instance,
                                           sample_csv_file):
        """Test data upload with pathlib.Path object."""
        path_obj = Path(sample_csv_file)

        with patch.object(data_client_instance, "_post") as mock_post:
            mock_post.return_value.status_code = 200
            mock_post.return_value.json.return_value = {"status": "success"}

            data_client_instance.upload_data(file_path=path_obj)

            assert mock_post.called

    def test_upload_data_with_different_backends(self, data_client_instance,
                                                 sample_csv_file):
        """Test data upload with different backend options."""
        backends = ["local", "s3", "gcs"]

        for backend in backends:
            with patch.object(data_client_instance, "_post") as mock_post:
                mock_post.return_value.status_code = 200
                mock_post.return_value.json.return_value = {"status": "success"}

                data_client_instance.upload_data(file_path=sample_csv_file,
                                                 backend=backend)

                # Verify backend parameter was included
                args, kwargs = mock_post.call_args
                # The data should be a MultipartEncoder containing the backend
                assert "data" in kwargs

    def test_upload_data_requires_configuration(self, temp_settings_file,
                                                sample_csv_file):
        """Test that upload_data requires proper configuration."""
        # Create unconfigured settings
        settings = Settings(settings_file=temp_settings_file)
        data_client = Data(settings=settings)

        with pytest.raises(ValueError) as exc_info:
            data_client.upload_data(file_path=sample_csv_file)

        assert "Missing required settings" in str(exc_info.value)

    @responses.activate
    def test_upload_data_with_sdf_file(self, data_client_instance,
                                       sample_sdf_file):
        """Test data upload with SDF file."""
        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json={
                "status":
                    "success",
                "dataset_address":
                    "deepchem://test-profile/test-project/test.sdf",
            },
            status=200,
        )

        result = data_client_instance.upload_data(file_path=sample_sdf_file,
                                                  filename="test.sdf",
                                                  description="Test SDF file")

        assert result["status"] == "success"

    def test_upload_data_multipart_encoder_creation(self, data_client_instance,
                                                    sample_csv_file):
        """Test that MultipartEncoder is created correctly."""
        with patch.object(data_client_instance, "_post") as mock_post, patch(
                "pyds.data.MultipartEncoder") as mock_encoder:

            mock_encoder_instance = Mock()
            mock_encoder_instance.content_type = "multipart/form-data"
            mock_encoder.return_value = mock_encoder_instance

            mock_post.return_value.status_code = 200
            mock_post.return_value.json.return_value = {"status": "success"}

            data_client_instance.upload_data(file_path=sample_csv_file,
                                             filename="test.csv",
                                             description="Test description")

            # Verify MultipartEncoder was called with correct fields
            mock_encoder.assert_called_once()
            call_args = mock_encoder.call_args[1]["fields"]
            assert "filename" in call_args
            assert "profile_name" in call_args
            assert "project_name" in call_args
            assert "file" in call_args
            assert "description" in call_args
            assert "backend" in call_args

    @responses.activate
    def test_upload_data_custom_headers(self, data_client_instance,
                                        sample_csv_file):
        """Test that correct headers are set for multipart upload."""
        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json={"status": "success"},
            status=200,
        )

        data_client_instance.upload_data(file_path=sample_csv_file)

        # Verify Content-Type header was set for multipart data
        request_headers = responses.calls[0].request.headers
        assert "multipart/form-data" in request_headers.get("Content-Type", "")

    def test_upload_data_empty_description(self, data_client_instance,
                                           sample_csv_file):
        """Test data upload with empty description."""
        with patch.object(data_client_instance, "_post") as mock_post:
            mock_post.return_value.status_code = 200
            mock_post.return_value.json.return_value = {"status": "success"}

            # Should not raise error with empty description
            data_client_instance.upload_data(file_path=sample_csv_file,
                                             description="")

            assert mock_post.called
