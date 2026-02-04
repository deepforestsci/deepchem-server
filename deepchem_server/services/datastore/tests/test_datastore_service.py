"""Tests for DatastoreService."""

import tempfile

import pytest

from deepchem_server.services.datastore.server.service import DatastoreService


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

    def test_upload_directory_files(self, service):
        """Test uploading multiple files as a directory."""
        profile = "test_profile"
        project = "test_project"
        key = "my_dataset"
        files = [
            ("file1.txt", b"content1"),
            ("subdir/file2.txt", b"content2"),
            ("subdir/nested/file3.txt", b"content3"),
        ]
        card = {"file_type": "dir", "data_type": "DiskDataset"}

        address = service.upload_directory_files(profile, project, key, files, card=card)
        assert address == f"deepchem://{profile}/{project}/{key}"

        # Verify files exist
        assert service.exists(profile, project, f"{key}/file1.txt")
        assert service.exists(profile, project, f"{key}/subdir/file2.txt")
        assert service.exists(profile, project, f"{key}/subdir/nested/file3.txt")

        # Verify card was stored
        retrieved_card = service.get_card(profile, project, key)
        assert retrieved_card == card

    def test_list_directory_contents(self, service):
        """Test listing files in a directory."""
        profile = "test_profile"
        project = "test_project"
        key = "list_test"
        files = [
            ("a.txt", b"a"),
            ("b.txt", b"b"),
            ("sub/c.txt", b"c"),
        ]

        service.upload_directory_files(profile, project, key, files)
        contents = service.list_directory_contents(profile, project, key)

        assert len(contents) == 3
        assert "a.txt" in contents
        assert "b.txt" in contents
        assert "sub/c.txt" in contents

    def test_list_directory_not_found(self, service):
        """Test listing nonexistent directory."""
        with pytest.raises(FileNotFoundError):
            service.list_directory_contents("no", "such", "dir")

    def test_download_directory_as_zip(self, service):
        """Test downloading a directory as ZIP."""
        import io
        import zipfile

        profile = "test_profile"
        project = "test_project"
        key = "zip_test"
        files = [
            ("file1.txt", b"hello"),
            ("nested/file2.txt", b"world"),
        ]

        service.upload_directory_files(profile, project, key, files)
        zip_data = service.download_directory_as_zip(profile, project, key)

        # Verify ZIP content
        zip_buffer = io.BytesIO(zip_data)
        with zipfile.ZipFile(zip_buffer, "r") as zf:
            names = zf.namelist()
            assert "file1.txt" in names
            assert "nested/file2.txt" in names
            assert zf.read("file1.txt") == b"hello"
            assert zf.read("nested/file2.txt") == b"world"

    def test_download_directory_not_found(self, service):
        """Test downloading nonexistent directory."""
        with pytest.raises(FileNotFoundError):
            service.download_directory_as_zip("no", "such", "dir")
