"""
Unit tests for Data class.
"""

import json
from pathlib import Path
from typing import Any, Dict
from unittest.mock import patch

import pandas as pd
import pytest
import requests
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

    def test_list_data_success(self, temp_test_file: str, data_client: Data) -> None:
        """Test list_data success."""

        with open(temp_test_file, "w") as f:
            f.write("test,data\n1,2")

        result = data_client.upload_data(file_path=temp_test_file, filename="test_data1.csv")
        dataset_address1 = result["dataset_address"]
        assert dataset_address1 is not None

        result = data_client.upload_data(file_path=temp_test_file, filename="test_data2.csv")
        dataset_address2 = result["dataset_address"]
        assert dataset_address2 is not None

        result = data_client.list_data()
        assert dataset_address1 in result
        assert dataset_address2 in result

    def test_to_datastore_key_strips_deepchem_uri(self, data_client: Data) -> None:
        """deepchem:// addresses normalize to the datastore key."""
        uri = "deepchem://test_profile/test_project/subdir/file.csv"
        assert data_client._to_datastore_key(uri) == "subdir/file.csv"

    def test_to_datastore_key_passthrough_plain_key(self, data_client: Data) -> None:
        """Bare filenames are returned unchanged."""
        assert data_client._to_datastore_key("metrics.json") == "metrics.json"

    @responses.activate
    def test_get_csv_parsed_as_dataframe(self, data_client: Data) -> None:
        """get parses CSV responses into a pandas DataFrame."""
        csv_body = b"a,b\n1,2\n3,4\n"
        responses.add(
            responses.GET,
            "http://localhost:8000/data/out.csv?profile_name=test_profile&project_name=test_project",
            body=csv_body,
            status=200,
            content_type="text/csv",
        )
        result = data_client.get("out.csv")
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ["a", "b"]
        assert len(result) == 2

    @responses.activate
    def test_get_normalizes_deepchem_address_in_request(self, data_client: Data) -> None:
        """get requests /data/<key> after stripping deepchem:// prefix."""
        responses.add(
            responses.GET,
            "http://localhost:8000/data/out.csv?profile_name=test_profile&project_name=test_project",
            body=b"x\n1\n",
            status=200,
        )
        addr = "deepchem://test_profile/test_project/out.csv"
        df = data_client.get(addr)
        assert isinstance(df, pd.DataFrame)
        assert list(df.columns) == ["x"]

    @responses.activate
    def test_get_json_parsed(self, data_client: Data) -> None:
        """get parses JSON into native Python objects."""
        payload = {"pearson_r2_score": 0.5, "rms_score": 0.1}
        responses.add(
            responses.GET,
            "http://localhost:8000/data/eval.json?profile_name=test_profile&project_name=test_project",
            body=json.dumps(payload).encode("utf-8"),
            status=200,
            content_type="application/json",
        )
        result = data_client.get("eval.json")
        assert result == payload

    @responses.activate
    def test_get_txt_as_string(self, data_client: Data) -> None:
        """get decodes plain-text extensions as str."""
        responses.add(
            responses.GET,
            "http://localhost:8000/data/readme.txt?profile_name=test_profile&project_name=test_project",
            body=b"hello\nworld",
            status=200,
        )
        assert data_client.get("readme.txt") == "hello\nworld"

    @responses.activate
    def test_get_writes_destination_path(self, data_client: Data, tmp_path: Path) -> None:
        """When destination_path is set, raw bytes are written to disk."""
        dest = tmp_path / "copy.csv"
        body = b"a\n1\n"
        responses.add(
            responses.GET,
            "http://localhost:8000/data/source.csv?profile_name=test_profile&project_name=test_project",
            body=body,
            status=200,
        )
        df = data_client.get("source.csv", destination_path=str(dest))
        assert dest.read_bytes() == body
        assert isinstance(df, pd.DataFrame)

    @responses.activate
    def test_get_http_error(self, data_client: Data) -> None:
        """Non-success status codes raise HTTPError."""
        responses.add(
            responses.GET,
            "http://localhost:8000/data/missing.csv?profile_name=test_profile&project_name=test_project",
            json={"detail": "Not Found"},
            status=404,
        )
        with pytest.raises(requests.exceptions.HTTPError):
            data_client.get("missing.csv")

    @responses.activate
    def test_get_unknown_extension_returns_message(self, data_client: Data) -> None:
        """Unsupported extensions yield the fallback string (no parser)."""
        responses.add(
            responses.GET,
            "http://localhost:8000/data/blob.bin?profile_name=test_profile&project_name=test_project",
            body=b"\x00\xff",
            status=200,
            content_type="application/octet-stream",
        )
        result = data_client.get("blob.bin")
        assert isinstance(result, str)
        assert "Unable to get data" in result

    @responses.activate
    def test_download_data_normalizes_deepchem_uri(self, data_client: Data, tmp_path: Path) -> None:
        """download_data uses the same key normalization as get."""
        dest = tmp_path / "dl.csv"
        body = b"x,y\n0,1\n"
        responses.add(
            responses.GET,
            "http://localhost:8000/data/out.csv?profile_name=test_profile&project_name=test_project",
            body=body,
            status=200,
        )
        uri = "deepchem://test_profile/test_project/out.csv"
        path = data_client.download_data(uri, str(dest))
        assert path == str(dest.resolve())
        assert dest.read_bytes() == body
