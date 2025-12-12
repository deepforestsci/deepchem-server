"""Tests for DatastoreService."""

import tempfile

import pytest

from deepchem_server.services.datastore.service import DatastoreService


@pytest.fixture
def temp_dir():
    """Create a temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def service(temp_dir):
    """Create a DatastoreService instance."""
    return DatastoreService(temp_dir)


class TestDatastoreService:
    """Tests for DatastoreService class."""

    def test_upload_and_get_data(self, service):
        """Test uploading and retrieving data."""
        profile = "test_profile"
        project = "test_project"
        key = "test_file.txt"
        data = b"Hello, World!"

        # Upload
        address = service.upload_data(profile, project, key, data)
        assert address == f"deepchem://{profile}/{project}/{key}"

        # Get
        retrieved = service.get_data(profile, project, key)
        assert retrieved == data

    def test_upload_with_card(self, service):
        """Test uploading with metadata card."""
        profile = "test_profile"
        project = "test_project"
        key = "data.csv"
        data = b"col1,col2\n1,2\n3,4"
        card = {"file_type": "csv", "rows": 2}

        service.upload_data(profile, project, key, data, card=card)

        retrieved_card = service.get_card(profile, project, key)
        assert retrieved_card == card

    def test_delete_object(self, service):
        """Test deleting an object."""
        profile = "test_profile"
        project = "test_project"
        key = "to_delete.txt"

        service.upload_data(profile, project, key, b"data")
        assert service.exists(profile, project, key)

        service.delete_object(profile, project, key)
        assert not service.exists(profile, project, key)

    def test_list_data(self, service):
        """Test listing data in a project."""
        profile = "test_profile"
        project = "test_project"

        service.upload_data(profile, project, "file1.txt", b"1")
        service.upload_data(profile, project, "file2.txt", b"2")
        service.upload_data(profile, project, "subdir/file3.txt", b"3")

        keys = service.list_data(profile, project)
        assert len(keys) == 3
        assert "file1.txt" in keys
        assert "file2.txt" in keys

    def test_exists(self, service):
        """Test exists check."""
        profile = "test_profile"
        project = "test_project"

        assert not service.exists(profile, project, "nonexistent.txt")

        service.upload_data(profile, project, "exists.txt", b"data")
        assert service.exists(profile, project, "exists.txt")

    def test_get_file_size(self, service):
        """Test getting file size."""
        profile = "test_profile"
        project = "test_project"
        key = "sized.txt"
        data = b"12345"

        service.upload_data(profile, project, key, data)
        size = service.get_file_size(profile, project, key)
        assert size == len(data)

    def test_get_nonexistent_raises(self, service):
        """Test that getting nonexistent file raises error."""
        with pytest.raises(FileNotFoundError):
            service.get_data("no", "such", "file.txt")

    def test_nested_keys(self, service):
        """Test keys with nested paths."""
        profile = "test_profile"
        project = "test_project"
        key = "path/to/nested/file.txt"
        data = b"nested content"

        address = service.upload_data(profile, project, key, data)
        assert address == f"deepchem://{profile}/{project}/{key}"

        retrieved = service.get_data(profile, project, key)
        assert retrieved == data
