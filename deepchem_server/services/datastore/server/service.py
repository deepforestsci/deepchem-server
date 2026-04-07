"""
DatastoreService - Core storage operations for DeepchemDatastore.

This service wraps file storage operations with thread-safe access and
follows the DeepchemAddress pattern (deepchem://profile/project/key).
"""

import json
import logging
import os
import shutil
import threading
from pathlib import Path
from typing import Any, Dict, List, Optional


logger = logging.getLogger(__name__)


class DatastoreService:
    """Thread-safe datastore service for file operations.

    Provides CRUD operations for data and model files following the
    DeepchemAddress structure: deepchem://profile/project/key

    Parameters
    ----------
    base_dir : str
        Base directory for file storage.
    """

    def __init__(self, base_dir: str) -> None:
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)
        self._lock = threading.RLock()
        logger.info(f"DatastoreService initialized with base_dir: {self.base_dir}")

    def _get_storage_path(self, profile: str, project: str, key: str) -> Path:
        """Get the full storage path for a given address."""
        return self.base_dir / profile / project / key

    def _get_card_path(self, profile: str, project: str, key: str, kind: str = "data") -> Path:
        """Get the card file path."""
        ext = ".cdc" if kind == "data" else ".cmc"
        return self.base_dir / profile / project / (key + ext)

    def upload_data(
        self,
        profile: str,
        project: str,
        key: str,
        data: bytes,
        card: Optional[Dict[str, Any]] = None,
        kind: str = "data",
    ) -> str:
        """Upload data to the datastore.

        Parameters
        ----------
        profile : str
            Profile name (e.g., company/user)
        project : str
            Project name
        key : str
            File key/name
        data : bytes
            File content as bytes
        card : dict, optional
            Metadata card for the file
        kind : str
            Type of object ('data' or 'model')

        Returns
        -------
        str
            DeepchemAddress of the uploaded file
        """
        with self._lock:
            storage_path = self._get_storage_path(profile, project, key)

            # Create parent directories
            storage_path.parent.mkdir(parents=True, exist_ok=True)

            # Write data
            storage_path.write_bytes(data)
            logger.info(f"Uploaded data to: {storage_path}")

            # Write card if provided
            if card is not None:
                card_path = self._get_card_path(profile, project, key, kind)
                card_path.write_text(json.dumps(card, indent=2))
                logger.debug(f"Wrote card to: {card_path}")

            return f"deepchem://{profile}/{project}/{key}"

    def create_directory(self, profile: str, project: str, dirname: str) -> str:
        """Create a directory in the datastore."""
        with self._lock:
            directory_path = self._get_storage_path(profile, project, dirname)
            directory_path.mkdir(parents=True, exist_ok=True)
            logger.info(f"Created directory: {directory_path}")
            return f"deepchem://{profile}/{project}/{dirname}"

    def move_object(self, profile: str, project: str, key: str, destination: str) -> str:
        """Move an object to a new location.

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name
        key : str
            File key/name
        destination : str
            New location address
        """
        with self._lock:
            storage_path = self._get_storage_path(profile, project, key)

            if not self.exists(profile, project, key):
                raise FileNotFoundError(f"Object not found: {storage_path}")

            prefix = "deepchem://"
            destination_profile = destination.replace(prefix, "").split("/")[0]
            destination_project = destination.replace(prefix, "").split("/")[1]
            destination_key = destination.replace(prefix + destination_profile + "/" + destination_project + "/", "")

            destination_path = self._get_storage_path(destination_profile, destination_project, destination_key)

            card = self.get_card(profile, project, key)
            if card:
                kind = "data" if "data_type" in card else "model" if "model_type" in card else "dir"
            else:
                kind = "dir"

            is_dir_object = (kind == "dir" or (kind == "data" and card is not None and card.get("file_type") == "dir"))
            if is_dir_object:
                self.create_directory(destination_profile, destination_project, destination_key)
                for file in self.list_directory_contents(profile, project, key):
                    self.move_object(
                        profile,
                        project,
                        f"deepchem://{profile}/{project}/{file}",
                        f"deepchem://{destination_profile}/{destination_project}/{file}",
                    )
            elif kind == "model":
                # move all the model files to the destination
                base_path = self._get_storage_path(profile, project, key)
                os.makedirs(destination_path, exist_ok=True)
                for file in os.listdir(base_path):
                    file_path = base_path / file
                    if file_path.is_file():
                        self.move_object(
                            profile,
                            project,
                            f"deepchem://{profile}/{project}/{file_path}",
                            f"deepchem://{destination_profile}/{destination_project}/{file_path}",
                        )
            else:
                self.upload_data(
                    destination_profile,
                    destination_project,
                    destination_key,
                    self.get_data(profile, project, key),
                )

            # move the card
            if card:
                card_path = self._get_card_path(profile, project, key, kind)
                card_path.rename(self._get_card_path(destination_profile, destination_project, destination_key, kind))

            self.delete_object(profile, project, key, kind)

            logger.info(f"Moved object from {storage_path} to {destination_path}")
            return destination

    def get_data(self, profile: str, project: str, key: str, fetch_sample: bool = False) -> bytes:
        """Get data from the datastore.

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name
        key : str
            File key/name
        fetch_sample : bool
            Whether to fetch only a sample (for large files)

        Returns
        -------
        bytes
            File content

        Raises
        ------
        FileNotFoundError
            If the file does not exist
        """
        storage_path = self._get_storage_path(profile, project, key)

        if not storage_path.exists():
            raise FileNotFoundError(f"Object not found: {storage_path}")

        if storage_path.is_dir():
            # For directories, return a JSON list of contents
            contents = [str(p.relative_to(storage_path)) for p in storage_path.rglob("*") if p.is_file()]
            return json.dumps(contents).encode()

        return storage_path.read_bytes()

    def delete_object(self, profile: str, project: str, key: str, kind: str = "data") -> bool:
        """Delete an object from the datastore.

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name
        key : str
            File key/name
        kind : str
            Type of object ('data', 'model', 'dir')

        Returns
        -------
        bool
            True if deletion was successful
        """
        with self._lock:
            storage_path = self._get_storage_path(profile, project, key)

            if storage_path.is_file():
                storage_path.unlink()
            elif storage_path.is_dir():
                shutil.rmtree(storage_path)
            else:
                raise FileNotFoundError(f"Object not found: {storage_path}")

            # Delete card if exists
            if kind != "dir":
                card_path = self._get_card_path(profile, project, key, kind)
                if card_path.exists():
                    card_path.unlink()

            logger.info(f"Deleted object: {storage_path}")
            return True

    def list_data(self, profile: str, project: str, include_card_files: bool = False) -> List[str]:
        """List all data keys in a profile/project (excludes card files).

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name
        include_card_files : bool
            Whether to include card files (.cdc, .cmc)
        Returns
        -------
        list of str
            List of keys in the project (without .cdc/.cmc files)
        """
        project_path = self.base_dir / profile / project

        if not project_path.exists():
            return []

        keys = []
        for path in project_path.rglob("*"):
            if path.is_file() and (include_card_files or path.suffix not in (".cdc", ".cmc")):
                keys.append(str(path.relative_to(project_path)))

        return sorted(keys)

    def list_all_objects(self, profile: str, project: str, prefix: str = "") -> List[str]:
        """List all objects including card files (.cdc, .cmc).

        This is used by DeepchemDatastore._get_datastore_objects() for
        checkpoint scanning during multicore featurization restart.

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name

        Returns
        -------
        list of str
            List of all file paths relative to project (including cards)
        """
        project_path = self.base_dir / profile / project

        if not project_path.exists():
            return []

        objects = []
        for path in project_path.rglob("*"):
            if path.is_file() and path.name.startswith(prefix):
                objects.append(str(path.relative_to(project_path)))
            elif path.is_dir():
                objects.append(str(path.relative_to(project_path)) + "/")
        return sorted(objects)

    def get_card(self, profile: str, project: str, key: str, kind: str = "data") -> Optional[Dict[str, Any]]:
        """Get the metadata card for an object.

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name
        key : str
            File key/name
        kind : str
            Type of object ('data' or 'model')

        Returns
        -------
        dict or None
            Card metadata, or None if not found
        """
        card_path = self._get_card_path(profile, project, key, kind)

        if not card_path.exists():
            return None

        return json.loads(card_path.read_text())

    def exists(self, profile: str, project: str, key: str) -> bool:
        """Check if an object exists.

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name
        key : str
            File key/name

        Returns
        -------
        bool
            True if the object exists
        """
        storage_path = self._get_storage_path(profile, project, key)
        return storage_path.exists()

    def get_file_size(self, profile: str, project: str, key: str) -> int:
        """Get the size of an object in bytes.

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name
        key : str
            File key/name

        Returns
        -------
        int
            Size in bytes
        """
        storage_path = self._get_storage_path(profile, project, key)

        if not storage_path.exists():
            raise FileNotFoundError(f"Object not found: {storage_path}")

        if storage_path.is_file():
            return storage_path.stat().st_size

        # For directories, sum all file sizes
        total_size = 0
        for path in storage_path.rglob("*"):
            if path.is_file():
                total_size += path.stat().st_size
        return total_size

    def upload_directory_files(
        self,
        profile: str,
        project: str,
        key: str,
        files: List[tuple],
        card: Optional[Dict[str, Any]] = None,
        kind: str = "data",
    ) -> str:
        """Upload multiple files as a directory structure.

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name
        key : str
            Directory key/name (base path for all files)
        files : list of tuples
            List of (relative_path, data_bytes) tuples
        card : dict, optional
            Metadata card for the directory
        kind : str
            Type of object ('data' or 'model')

        Returns
        -------
        str
            DeepchemAddress of the uploaded directory
        """
        with self._lock:
            base_path = self._get_storage_path(profile, project, key)
            base_path.mkdir(parents=True, exist_ok=True)

            for relative_path, data in files:
                file_path = base_path / relative_path
                file_path.parent.mkdir(parents=True, exist_ok=True)
                file_path.write_bytes(data)
                logger.debug(f"Uploaded file: {file_path}")

            logger.info(f"Uploaded directory with {len(files)} files to: {base_path}")

            # Write card if provided
            if card is not None:
                card_path = self._get_card_path(profile, project, key, kind)
                card_path.write_text(json.dumps(card, indent=2))
                logger.debug(f"Wrote card to: {card_path}")

            return f"deepchem://{profile}/{project}/{key}"

    def list_directory_contents(self, profile: str, project: str, key: str) -> List[str]:
        """List all files in a directory.

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name
        key : str
            Directory key/name

        Returns
        -------
        list of str
            List of relative file paths within the directory
        """
        storage_path = self._get_storage_path(profile, project, key)

        if not storage_path.exists():
            raise FileNotFoundError(f"Directory not found: {storage_path}")

        if not storage_path.is_dir():
            raise ValueError(f"Path is not a directory: {storage_path}")

        files = []
        for path in storage_path.rglob("*"):
            if path.is_file():
                files.append(str(path.relative_to(storage_path)))

        return sorted(files)

    def download_directory_as_zip(self, profile: str, project: str, key: str) -> bytes:
        """Download a directory as a ZIP archive.

        Creates the ZIP on-demand for user downloads.

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name
        key : str
            Directory key/name

        Returns
        -------
        bytes
            ZIP archive content
        """
        import io
        import zipfile

        storage_path = self._get_storage_path(profile, project, key)

        if not storage_path.exists():
            raise FileNotFoundError(f"Directory not found: {storage_path}")

        if not storage_path.is_dir():
            raise ValueError(f"Path is not a directory: {storage_path}")

        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
            for file_path in storage_path.rglob("*"):
                if file_path.is_file():
                    arcname = str(file_path.relative_to(storage_path))
                    zf.write(file_path, arcname)

        logger.info(f"Created ZIP archive for directory: {storage_path}")
        return zip_buffer.getvalue()
