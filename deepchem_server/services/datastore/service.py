"""
DatastoreService - Core storage operations for DeepchemDatastore.

This service wraps file storage operations with thread-safe access and
follows the DeepchemAddress pattern (deepchem://profile/project/key).
"""

import json
import logging
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
            contents = [
                str(p.relative_to(storage_path)) for p in storage_path.rglob("*") if p.is_file()
            ]
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

    def list_data(self, profile: str, project: str) -> List[str]:
        """List all data keys in a profile/project.

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name

        Returns
        -------
        list of str
            List of keys in the project
        """
        project_path = self.base_dir / profile / project

        if not project_path.exists():
            return []

        keys = []
        for path in project_path.rglob("*"):
            if path.is_file() and not path.suffix in (".cdc", ".cmc"):
                keys.append(str(path.relative_to(project_path)))

        return sorted(keys)

    def get_card(
        self, profile: str, project: str, key: str, kind: str = "data"
    ) -> Optional[Dict[str, Any]]:
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
