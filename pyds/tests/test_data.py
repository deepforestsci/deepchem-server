"""
Unit tests for Data class.
"""

from pathlib import Path
from typing import Any, Dict
from unittest.mock import patch

import pytest
import responses

from pyds.data import Data
from pyds.settings import Settings


class TestData:
    """Unit tests for Data class."""

    def test_init(self, test_settings: Settings) -> None:
        """Test Data client initialization."""
        client = Data(settings=test_settings)

        assert client.settings == test_settings
        assert client.base_url == "http://localhost:8000"

    def test_upload_data_file_not_found(self, data_client: Data) -> None:
        """Test upload_data with non-existent file."""
        with pytest.raises(FileNotFoundError, match="File not found"):
            data_client.upload_data("nonexistent_file.csv")

    @responses.activate
    def test_upload_data_success(self, data_client: Data, temp_test_file: str,
                                 sample_upload_response: Dict[str, Any]) -> None:
        """Test successful data upload."""
        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json=sample_upload_response,
            status=200,
        )

        with open(temp_test_file, "w") as f:
            f.write("test,data\n1,2")

        result = data_client.upload_data(file_path=temp_test_file,
                                         filename="custom_name.csv",
                                         description="Test description")

        assert result == sample_upload_response

    @responses.activate
    def test_upload_data_with_defaults(self, data_client: Data, temp_test_file: str,
                                       sample_upload_response: Dict[str, Any]) -> None:
        """Test data upload with default parameters."""
        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json=sample_upload_response,
            status=200,
        )

        with open(temp_test_file, "w") as f:
            f.write("test,data\n1,2")

        result = data_client.upload_data(file_path=temp_test_file)

        assert result == sample_upload_response

    @responses.activate
    def test_upload_data_with_profile_project_override(self, data_client: Data, temp_test_file: str,
                                                       sample_upload_response: Dict[str, Any]) -> None:
        """Test data upload with profile and project override."""
        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json=sample_upload_response,
            status=200,
        )

        # Write actual content to temp file
        with open(temp_test_file, "w") as f:
            f.write("test,data\n1,2")

        result = data_client.upload_data(
            file_path=temp_test_file,
            profile_name="custom_profile",
            project_name="custom_project",
        )

        assert result == sample_upload_response

    def test_upload_data_missing_settings(self, test_settings_not_configured: Settings, temp_test_file: str) -> None:
        """Test upload_data with missing settings."""
        client = Data(settings=test_settings_not_configured)

        with pytest.raises(ValueError, match="Missing required settings"):
            client.upload_data(file_path=temp_test_file)

    @responses.activate
    def test_upload_data_api_error(self, data_client: Data, temp_test_file: str) -> None:
        """Test upload_data with API error."""
        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json={"detail": "Upload failed"},
            status=400,
        )

        with open(temp_test_file, "w") as f:
            f.write("test,data\n1,2")

        with pytest.raises(ValueError, match="Failed to upload data"):
            data_client.upload_data(file_path=temp_test_file)

    def test_upload_data_file_path_as_string(self, data_client: Data, temp_test_file: str) -> None:
        """Test upload_data handles string file path."""
        file_path_str = str(temp_test_file)

        with patch.object(Path, "exists", return_value=False):
            with pytest.raises(FileNotFoundError):
                data_client.upload_data(file_path=file_path_str)

    @responses.activate
    def test_upload_data_filename_from_path(self, data_client: Data, temp_test_file: str) -> None:
        """Test upload_data uses filename from path when not provided."""

        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json={"dataset_address": "test"},
            status=200,
        )

        with open(temp_test_file, "w") as f:
            f.write("test,data\n1,2")

        result = data_client.upload_data(file_path=temp_test_file)

        assert result == {"dataset_address": "test"}

    @responses.activate
    def test_upload_data_request_exception(self, data_client: Data, temp_test_file: str) -> None:
        """Test upload_data with request exception."""
        with open(temp_test_file, "w") as f:
            f.write("test,data\n1,2")

        with patch.object(data_client, "_post", side_effect=Exception("Network error")):
            with pytest.raises(ValueError, match="Failed to upload data: Network error"):
                data_client.upload_data(file_path=temp_test_file)

    @responses.activate
    def test_upload_data_file_closes_on_success(self, data_client: Data, temp_test_file: str) -> None:
        """Test that file handle is closed on successful upload."""
        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json={"dataset_address": "test"},
            status=200,
        )

        with open(temp_test_file, "w") as f:
            f.write("test,data\n1,2")

        result = data_client.upload_data(file_path=temp_test_file)
        assert result == {"dataset_address": "test"}

    def test_upload_data_file_closes_on_error(self, data_client: Data, temp_test_file: str) -> None:
        """Test that file handle is closed on upload error."""
        with open(temp_test_file, "w") as f:
            f.write("test,data\n1,2")

        with patch.object(data_client, "_post", side_effect=Exception("Upload error")):
            with pytest.raises(ValueError):
                data_client.upload_data(file_path=temp_test_file)

    @responses.activate
    def test_upload_data_multipart_encoder_fields(self, data_client: Data, temp_test_file: str) -> None:
        """Test upload_data with all parameters."""
        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json={"dataset_address": "test"},
            status=200,
        )

        with open(temp_test_file, "w") as f:
            f.write("test,data\n1,2")

        result = data_client.upload_data(
            file_path=temp_test_file,
            filename="test.csv",
            description="Test file",
            backend="custom",
        )

        assert result == {"dataset_address": "test"}

    @responses.activate
    def test_upload_data_no_description(self, data_client: Data, temp_test_file: str) -> None:
        """Test upload_data without description."""
        responses.add(
            responses.POST,
            "http://localhost:8000/data/uploaddata",
            json={"dataset_address": "test"},
            status=200,
        )

        with open(temp_test_file, "w") as f:
            f.write("test,data\n1,2")

        result = data_client.upload_data(file_path=temp_test_file)

        assert result == {"dataset_address": "test"}
