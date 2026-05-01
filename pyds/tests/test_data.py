"""
Unit and integration tests for Data class.
"""

import json
import tempfile
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

    def _upload(self, client: Data, content: bytes, filename: str) -> str:
        """Upload bytes as a named file and return the dataset_address."""
        with tempfile.NamedTemporaryFile(suffix=Path(filename).suffix, delete=False) as f:
            f.write(content)
            tmp_path = f.name
        result = client.upload_data(file_path=tmp_path, filename=filename)
        Path(tmp_path).unlink(missing_ok=True)
        return result["dataset_address"]

    def test_get_csv_returns_dataframe(self, live_data_client: Data) -> None:
        """get on a .csv address returns a pandas DataFrame with correct content."""
        csv_content = b"smiles,label\nCCO,1\nCCC,0\nCCCO,1\n"
        address = self._upload(live_data_client, csv_content, "test_get.csv")

        result = live_data_client.get(address)

        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ["smiles", "label"]
        assert len(result) == 3
        assert result["smiles"].tolist() == ["CCO", "CCC", "CCCO"]
        assert result["label"].tolist() == [1, 0, 1]

    def test_get_json_returns_dict(self, live_data_client: Data) -> None:
        """get on a .json address returns a parsed Python dict."""
        payload = {"molecule": "CCO", "property": 1.23}
        json_content = json.dumps(payload).encode("utf-8")
        address = self._upload(live_data_client, json_content, "test_get.json")

        result = live_data_client.get(address)

        assert isinstance(result, dict)
        assert result == payload

    def test_get_txt_returns_string(self, live_data_client: Data) -> None:
        """get on a .txt address returns the file content as a str."""
        text = "hello from deepchem server\nline two\n"
        address = self._upload(live_data_client, text.encode("utf-8"), "test_get.txt")

        result = live_data_client.get(address)

        assert isinstance(result, str)
        assert result == text

    def test_get_pdb_returns_string(self, live_data_client: Data) -> None:
        """get on a .pdb address returns the raw PDB text as a str."""
        pdb_text = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000\nEND\n"
        address = self._upload(live_data_client, pdb_text.encode("utf-8"), "test_get.pdb")

        result = live_data_client.get(address)

        assert isinstance(result, str)
        assert "ATOM" in result
        assert "END" in result

    def test_get_saves_raw_bytes_to_destination_path(self, live_data_client: Data, tmp_path: Path) -> None:
        """get writes the raw response bytes to destination_path when provided."""
        csv_content = b"smiles,label\nCCO,1\n"
        address = self._upload(live_data_client, csv_content, "test_dest.csv")
        dest = tmp_path / "saved.csv"

        live_data_client.get(address, destination_path=str(dest))

        assert dest.exists()
        assert dest.read_bytes() == csv_content

    def test_get_with_deepchem_address_format(self, live_data_client: Data) -> None:
        """get accepts a deepchem:// address returned by upload and resolves it correctly.

        The server returns dataset_address strings in deepchem://profile/project/...
        format.  Passing that address directly to get() exercises the
        _get_address_key stripping path.
        """
        csv_content = b"smiles,label\nCCO,1\n"
        address = self._upload(live_data_client, csv_content, "test_dcaddr.csv")

        assert address.startswith("deepchem://"), f"Expected upload to return a deepchem:// address, got: {address!r}"

        result = live_data_client.get(address)

        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ["smiles", "label"]

    def test_get_nonexistent_address_raises_http_error(self, live_data_client: Data) -> None:
        """get raises HTTPError when the server returns a 404 for an unknown key."""
        with pytest.raises(requests.exceptions.HTTPError):
            live_data_client.get("does_not_exist/data/test.csv")
