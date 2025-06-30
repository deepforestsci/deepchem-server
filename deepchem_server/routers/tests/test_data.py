import pytest
from unittest.mock import Mock, patch, AsyncMock
from fastapi.testclient import TestClient
from fastapi import UploadFile
from deepchem_server.main import app
from deepchem_server.core.cards import DataCard
from deepchem_server.core.datastore import DiskDataStore
import io


@pytest.fixture
def client():
    """Create a test client for the FastAPI app"""
    return TestClient(app)


@pytest.fixture
def mock_datastore():
    """Create a mock datastore for testing"""
    mock_store = Mock(spec=DiskDataStore)
    mock_store._get_datastore_objects.return_value = [
        "dataset1.csv",
        "dataset2.sdf",
        "experiment_data.txt",
    ]
    mock_store.storage_loc = "/mock/storage"
    return mock_store


class TestUploadData:
    """Test cases for POST /data/uploaddata endpoint"""

    @patch("deepchem_server.routers.data._upload_data")
    def test_upload_data_success_csv(self, mock_upload, client):
        """Test successful CSV file upload"""
        mock_upload.return_value = "deepchem://testprofile/testproject/test.csv"

        # Create test file content
        file_content = b"smiles,label\nCCC,1\nCCCC,0"

        response = client.post(
            "/data/uploaddata",
            data={
                "profile_name": "testprofile",
                "project_name": "testproject",
                "filename": "test.csv",
                "description": "Test CSV file",
            },
            files={"file": ("test.csv", io.BytesIO(file_content), "text/csv")},
        )

        assert response.status_code == 200
        assert response.json() == {"dataset_address": "deepchem://testprofile/testproject/test.csv"}

        # Verify mock was called with correct parameters
        mock_upload.assert_called_once()
        call_args = mock_upload.call_args
        assert call_args[0][0] == "testprofile"  # profile_name
        assert call_args[0][1] == "testproject"  # project_name
        assert call_args[0][2] == "test.csv"  # filename
        assert call_args[0][3] == file_content  # contents

        # Verify DataCard was created properly
        data_card = call_args[0][4]
        assert isinstance(data_card, DataCard)
        assert data_card.file_type == "csv"
        assert data_card.data_type == "pandas.DataFrame"
        assert data_card.description == "Test CSV file"

    @patch("deepchem_server.routers.data._upload_data")
    def test_upload_data_no_filename_uses_original(self, mock_upload, client):
        """Test upload without filename uses original file.filename"""
        mock_upload.return_value = "deepchem://testprofile/testproject/original.sdf"

        file_content = b"mol file content"

        response = client.post(
            "/data/uploaddata",
            data={
                "profile_name": "testprofile",
                "project_name": "testproject",
                "description": "Test SDF file",
            },
            files={"file": ("original.sdf", io.BytesIO(file_content), "chemical/x-mdl-sdfile")},
        )

        assert response.status_code == 200

        # Verify filename was taken from uploaded file
        call_args = mock_upload.call_args
        assert call_args[0][2] == "original.sdf"

    @patch("deepchem_server.routers.data._upload_data")
    def test_upload_data_different_file_types(self, mock_upload, client):
        """Test upload with different file types creates correct data cards"""
        test_cases = [
            ("test.parquet", "pandas.DataFrame"),
            ("test.pdb", "text/plain"),
            ("test.zip", "binary"),
            ("test.png", "png"),
            ("test.unknown", ""),
        ]

        for filename, expected_data_type in test_cases:
            mock_upload.return_value = f"deepchem://test/test/{filename}"

            response = client.post(
                "/data/uploaddata",
                data={
                    "profile_name": "test",
                    "project_name": "test",
                    "filename": filename,
                    "description": "Test file",
                },
                files={"file": (filename, io.BytesIO(b"content"), "application/octet-stream")},
            )

            assert response.status_code == 200

            # Check DataCard type inference
            call_args = mock_upload.call_args
            data_card = call_args[0][4]
            assert data_card.data_type == expected_data_type

    @patch("deepchem_server.routers.data._upload_data")
    def test_upload_data_dict_description(self, mock_upload, client):
        """Test upload with dictionary description (passed through as-is)"""
        mock_upload.return_value = "deepchem://test/test/test.csv"

        # When description is not a string, it should be used as-is for the card
        description_dict = {"custom": "metadata", "version": "1.0"}

        response = client.post(
            "/data/uploaddata",
            data={
                "profile_name": "test",
                "project_name": "test",
                "filename": "test.csv",
                "description": str(description_dict),  # FastAPI form data converts to string
            },
            files={"file": ("test.csv", io.BytesIO(b"data"), "text/csv")},
        )

        assert response.status_code == 200

    @patch("deepchem_server.routers.data._upload_data")
    def test_upload_data_failure(self, mock_upload, client):
        """Test upload failure when _upload_data returns None"""
        mock_upload.return_value = None

        response = client.post(
            "/data/uploaddata",
            data={
                "profile_name": "test",
                "project_name": "test",
                "filename": "test.csv",
                "description": "Test file",
            },
            files={"file": ("test.csv", io.BytesIO(b"data"), "text/csv")},
        )

        assert response.status_code == 500
        assert response.json() == {"detail": "Failed to upload data"}

    def test_upload_data_missing_file(self, client):
        """Test upload without file parameter"""
        response = client.post(
            "/data/uploaddata",
            data={"profile_name": "test", "project_name": "test", "description": "Test file"},
        )

        assert response.status_code == 422  # Unprocessable Entity due to missing required field


class TestListDatastore:
    """Test cases for GET /data/list endpoint"""

    @patch("deepchem_server.routers.data._init_datastore")
    def test_list_datastore_success(self, mock_init_datastore, client, mock_datastore):
        """Test successful datastore listing"""
        mock_init_datastore.return_value = mock_datastore

        response = client.get(
            "/data/list", params={"profile_name": "testprofile", "project_name": "testproject"}
        )

        assert response.status_code == 200
        data = response.json()

        assert data["count"] == 3
        assert data["profile"] == "testprofile"
        assert data["project"] == "testproject"
        assert len(data["objects"]) == 3
        assert "deepchem://testprofile/testproject/dataset1.csv" in data["objects"]
        assert "deepchem://testprofile/testproject/dataset2.sdf" in data["objects"]
        assert "deepchem://testprofile/testproject/experiment_data.txt" in data["objects"]

    @patch("deepchem_server.routers.data._init_datastore")
    def test_list_datastore_empty(self, mock_init_datastore, client):
        """Test listing empty datastore"""
        mock_datastore = Mock(spec=DiskDataStore)
        mock_datastore._get_datastore_objects.return_value = []
        mock_init_datastore.return_value = mock_datastore

        response = client.get(
            "/data/list", params={"profile_name": "testprofile", "project_name": "testproject"}
        )

        assert response.status_code == 200
        data = response.json()
        assert data["count"] == 0
        assert data["objects"] == []

    @patch("deepchem_server.routers.data._init_datastore")
    def test_list_datastore_unsupported_type(self, mock_init_datastore, client):
        """Test error when datastore is not DiskDataStore"""
        mock_init_datastore.return_value = Mock()  # Not a DiskDataStore

        response = client.get(
            "/data/list", params={"profile_name": "testprofile", "project_name": "testproject"}
        )

        assert response.status_code == 500
        assert "Unsupported datastore type" in response.json()["detail"]

    @patch("deepchem_server.routers.data._init_datastore")
    def test_list_datastore_exception(self, mock_init_datastore, client):
        """Test error handling when datastore operations fail"""
        mock_init_datastore.side_effect = Exception("Datastore connection failed")

        response = client.get(
            "/data/list", params={"profile_name": "testprofile", "project_name": "testproject"}
        )

        assert response.status_code == 500
        assert "Failed to list datastore" in response.json()["detail"]

    def test_list_datastore_missing_params(self, client):
        """Test error when required parameters are missing"""
        response = client.get("/data/list")

        assert response.status_code == 422  # Missing required query parameters

    @patch("deepchem_server.routers.data._init_datastore")
    def test_list_datastore_custom_backend(self, mock_init_datastore, client, mock_datastore):
        """Test listing with custom backend parameter"""
        mock_init_datastore.return_value = mock_datastore

        response = client.get(
            "/data/list",
            params={
                "profile_name": "testprofile",
                "project_name": "testproject",
                "backend": "custom",
            },
        )

        assert response.status_code == 200
        # Verify backend parameter was passed to _init_datastore
        mock_init_datastore.assert_called_once_with(
            profile_name="testprofile", project_name="testproject", backend="custom"
        )


class TestSearchDatastore:
    """Test cases for GET /data/search endpoint"""

    @patch("deepchem_server.routers.data._init_datastore")
    @patch("deepchem_server.routers.data.process.extract")
    def test_search_datastore_success(
        self, mock_extract, mock_init_datastore, client, mock_datastore
    ):
        """Test successful datastore search"""
        mock_init_datastore.return_value = mock_datastore
        # Mock fuzzy matching results: (match, score)
        mock_extract.return_value = [
            ("dataset1.csv", 85),
            ("experiment_data.txt", 70),
            ("dataset2.sdf", 45),  # Below default cutoff of 60
        ]

        response = client.get(
            "/data/search",
            params={
                "query": "dataset",
                "profile_name": "testprofile",
                "project_name": "testproject",
            },
        )

        assert response.status_code == 200
        data = response.json()

        assert data["count"] == 2  # Only results above score_cutoff=60
        assert data["query"] == "dataset"
        assert data["profile"] == "testprofile"
        assert data["project"] == "testproject"
        assert data["score_cutoff"] == 60

        results = data["results"]
        assert len(results) == 2

        # Check first result
        assert results[0]["object_name"] == "dataset1.csv"
        assert results[0]["similarity_score"] == 85
        assert results[0]["full_address"] == "deepchem://testprofile/testproject/dataset1.csv"

        # Check second result
        assert results[1]["object_name"] == "experiment_data.txt"
        assert results[1]["similarity_score"] == 70

    @patch("deepchem_server.routers.data._init_datastore")
    @patch("deepchem_server.routers.data.process.extract")
    def test_search_datastore_custom_params(
        self, mock_extract, mock_init_datastore, client, mock_datastore
    ):
        """Test search with custom limit and score_cutoff"""
        mock_init_datastore.return_value = mock_datastore
        mock_extract.return_value = [
            ("dataset1.csv", 90),
            ("dataset2.sdf", 80),
            ("experiment_data.txt", 75),
        ]

        response = client.get(
            "/data/search",
            params={
                "query": "data",
                "profile_name": "testprofile",
                "project_name": "testproject",
                "limit": 2,
                "score_cutoff": 80,
            },
        )

        assert response.status_code == 200
        data = response.json()

        assert data["count"] == 2  # Results above score_cutoff=80
        assert data["score_cutoff"] == 80

        # Verify extract was called with correct limit
        mock_extract.assert_called_once_with(
            "data", ["dataset1.csv", "dataset2.sdf", "experiment_data.txt"], limit=2
        )

    @patch("deepchem_server.routers.data._init_datastore")
    def test_search_datastore_empty(self, mock_init_datastore, client):
        """Test search on empty datastore"""
        mock_datastore = Mock(spec=DiskDataStore)
        mock_datastore._get_datastore_objects.return_value = []
        mock_init_datastore.return_value = mock_datastore

        response = client.get(
            "/data/search",
            params={
                "query": "anything",
                "profile_name": "testprofile",
                "project_name": "testproject",
            },
        )

        assert response.status_code == 200
        data = response.json()
        assert data["count"] == 0
        assert data["results"] == []

    @patch("deepchem_server.routers.data._init_datastore")
    @patch("deepchem_server.routers.data.process.extract")
    def test_search_datastore_no_matches_above_cutoff(
        self, mock_extract, mock_init_datastore, client, mock_datastore
    ):
        """Test search where no results meet score cutoff"""
        mock_init_datastore.return_value = mock_datastore
        mock_extract.return_value = [("dataset1.csv", 50), ("dataset2.sdf", 45)]

        response = client.get(
            "/data/search",
            params={
                "query": "nonexistent",
                "profile_name": "testprofile",
                "project_name": "testproject",
                "score_cutoff": 60,
            },
        )

        assert response.status_code == 200
        data = response.json()
        assert data["count"] == 0
        assert data["results"] == []

    @patch("deepchem_server.routers.data._init_datastore")
    def test_search_datastore_unsupported_type(self, mock_init_datastore, client):
        """Test error when datastore is not DiskDataStore"""
        mock_init_datastore.return_value = Mock()  # Not a DiskDataStore

        response = client.get(
            "/data/search",
            params={"query": "test", "profile_name": "testprofile", "project_name": "testproject"},
        )

        assert response.status_code == 500
        assert "Unsupported datastore type" in response.json()["detail"]

    @patch("deepchem_server.routers.data._init_datastore")
    def test_search_datastore_exception(self, mock_init_datastore, client):
        """Test error handling when search operations fail"""
        mock_init_datastore.side_effect = Exception("Search failed")

        response = client.get(
            "/data/search",
            params={"query": "test", "profile_name": "testprofile", "project_name": "testproject"},
        )

        assert response.status_code == 500
        assert "Failed to search datastore" in response.json()["detail"]

    def test_search_datastore_missing_query(self, client):
        """Test error when query parameter is missing"""
        response = client.get(
            "/data/search", params={"profile_name": "testprofile", "project_name": "testproject"}
        )

        assert response.status_code == 422  # Missing required query parameter
