"""
Python HTTP client for DeepchemDatastore service.

Provides a low-level interface for other services to communicate
with the datastore service via HTTP.
"""
import json
import logging
import os
from typing import Any, Dict, List, Optional
from urllib.parse import urlparse

import httpx


logger = logging.getLogger(__name__)


class DatastoreClient:
    """Low-level HTTP client for communicating with DeepchemDatastore service.

    Parameters
    ----------
    url : str, optional
        Full URL of the datastore service (e.g., "http://localhost:8081")
    host : str, optional
        Hostname of the datastore service (deprecated, use url)
    port : int, optional
        Port of the datastore service (deprecated, use url)
    api_key : str
        API key for authentication
    timeout : float
        Request timeout in seconds
    """

    def __init__(
        self,
        url: Optional[str] = None,
        host: str = "localhost",
        port: int = 8081,
        api_key: Optional[str] = None,
        timeout: float = 300.0,
    ) -> None:
        if url:
            # Parse URL and extract base
            parsed = urlparse(url)
            base = f"{parsed.scheme}://{parsed.netloc}"
            self.base_url = f"{base}/api/v1"
        else:
            self.base_url = f"http://{host}:{port}/api/v1"

        self.api_key = api_key or os.getenv("DATASTORE_API_KEY", "dev-api-key")
        self.timeout = timeout
        self._client: Optional[httpx.Client] = None

    def _get_client(self) -> httpx.Client:
        """Get or create the HTTP client."""
        if self._client is None:
            self._client = httpx.Client(
                timeout=self.timeout,
                headers={"X-API-Key": self.api_key},  # type: ignore
            )
        return self._client

    def close(self) -> None:
        """Close the HTTP client."""
        if self._client is not None:
            self._client.close()
            self._client = None

    def __enter__(self) -> "DatastoreClient":
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    def _parse_address(self, address: str) -> tuple:
        """Parse a deepchem address into profile, project, key."""
        if address.startswith("deepchem://"):
            address = address[len("deepchem://"):]
        parts = address.split("/", 2)
        if len(parts) < 3:
            raise ValueError(f"Invalid address format: {address}")
        return parts[0], parts[1], parts[2]

    def _make_address(self, profile: str, project: str, key: str) -> str:
        """Create a deepchem address from components."""
        return f"deepchem://{profile}/{project}/{key}"

    def upload_data(
        self,
        address: str,
        data: bytes,
        card: Optional[Dict[str, Any]] = None,
        kind: str = "data",
    ) -> str:
        """Upload data to the datastore.

        Parameters
        ----------
        address : str
            DeepchemAddress (deepchem://profile/project/key)
        data : bytes
            File content
        card : dict, optional
            Metadata card
        kind : str
            Type of object

        Returns
        -------
        str
            The assigned address
        """
        profile, project, key = self._parse_address(address)

        files = {"file": ("data", data)}
        form_data: Dict[str, Any] = {"kind": kind}
        if card:
            form_data["card"] = json.dumps(card)

        response = self._get_client().post(
            f"{self.base_url}/data/{profile}/{project}/{key}",
            files=files,
            data=form_data,
        )
        response.raise_for_status()
        return response.json()["address"]

    def create_directory(
        self,
        address: str,
    ) -> str:
        """Create a directory in the datastore."""
        profile, project, key = self._parse_address(address)
        response = self._get_client().post(f"{self.base_url}/dir/{profile}/{project}/{key}")
        response.raise_for_status()
        return response.json()["address"]

    def move_object(
        self,
        address: str,
        destination: str,
    ) -> str:
        """Move an object to a new location.

        Parameters
        ----------
        address : str
            DeepchemAddress
        destination : str
            New location address

        Returns
        -------
        str
            The assigned address
        """
        profile, project, key = self._parse_address(address)
        response = self._get_client().post(
            f"{self.base_url}/move/{profile}/{project}/{key}", json={"destination": destination}
        )
        response.raise_for_status()
        return response.json()["address"]

    def get_data(
        self,
        address: str,
        fetch_sample: bool = False,
    ) -> bytes:
        """Get data from the datastore.

        Parameters
        ----------
        address : str
            DeepchemAddress
        fetch_sample : bool
            Whether to fetch only a sample

        Returns
        -------
        bytes
            File content
        """
        profile, project, key = self._parse_address(address)

        response = self._get_client().get(
            f"{self.base_url}/data/{profile}/{project}/{key}",
            params={"fetch_sample": fetch_sample},
        )
        response.raise_for_status()
        return response.content

    def delete_object(
        self,
        address: str,
        kind: str = "data",
    ) -> bool:
        """Delete an object from the datastore."""
        profile, project, key = self._parse_address(address)

        response = self._get_client().delete(
            f"{self.base_url}/data/{profile}/{project}/{key}",
            params={"kind": kind},
        )
        response.raise_for_status()
        return response.json()["status"] == "success"

    def list_data(self, profile: str, project: str, include_card_files: bool = False) -> List[str]:
        """List all keys in a profile/project (excludes card files).

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name
        include_card_files : bool
            Whether to include card files (.cdc, .cmc)
        """
        response = self._get_client().get(
            f"{self.base_url}/list/{profile}/{project}",
            params={"include_card_files": include_card_files},
        )
        response.raise_for_status()
        return response.json()["keys"]

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
        prefix : str
            Prefix to filter objects
        Returns
        -------
        list of str
            All objects including .cdc and .cmc files
        """
        response = self._get_client().get(f"{self.base_url}/list-all/{profile}/{project}", params={"prefix": prefix})
        response.raise_for_status()
        return response.json()["objects"]

    def get_card(
        self,
        address: str,
        kind: str = "data",
    ) -> Optional[Dict[str, Any]]:
        """Get the metadata card for an object."""
        profile, project, key = self._parse_address(address)

        response = self._get_client().get(
            f"{self.base_url}/card/{profile}/{project}/{key}",
            params={"kind": kind},
        )
        response.raise_for_status()
        return response.json()["card"]

    def exists(self, address: str) -> bool:
        """Check if an object exists."""
        profile, project, key = self._parse_address(address)

        response = self._get_client().get(f"{self.base_url}/exists/{profile}/{project}/{key}")
        response.raise_for_status()
        return response.json()["exists"]

    def get_size(self, address: str) -> int:
        """Get the size of an object."""
        profile, project, key = self._parse_address(address)

        response = self._get_client().get(f"{self.base_url}/size/{profile}/{project}/{key}")
        response.raise_for_status()
        return response.json()["size"]

    def healthcheck(self) -> bool:
        """Check if the service is healthy."""
        try:
            response = httpx.get(f"{self.base_url}/healthcheck", timeout=5.0)
            return response.status_code == 200
        except Exception:
            return False

    def upload_file(
        self,
        address: str,
        file_path: str,
        card: Optional[Dict[str, Any]] = None,
        kind: str = "data",
    ) -> str:
        """Upload a file from disk to the datastore.

        Parameters
        ----------
        address : str
            DeepchemAddress
        file_path : str
            Path to file on disk
        card : dict, optional
            Metadata card
        kind : str
            Type of object

        Returns
        -------
        str
            The assigned address
        """
        with open(file_path, "rb") as f:
            data = f.read()
        return self.upload_data(address, data, card, kind)

    def download_file(
        self,
        address: str,
        dest_path: str,
        fetch_sample: bool = False,
    ) -> str:
        """Download data to a file on disk.

        Parameters
        ----------
        address : str
            DeepchemAddress
        dest_path : str
            Destination path on disk
        fetch_sample : bool
            Whether to fetch only a sample

        Returns
        -------
        str
            Path to downloaded file
        """
        data = self.get_data(address, fetch_sample)

        # Ensure parent directory exists
        os.makedirs(os.path.dirname(dest_path) or ".", exist_ok=True)

        with open(dest_path, "wb") as f:
            f.write(data)
        return dest_path

    def upload_directory(
        self,
        address: str,
        dir_path: str,
        card: Optional[Dict[str, Any]] = None,
        kind: str = "data",
    ) -> str:
        """Upload a directory by uploading each file directly.

        Files are uploaded individually preserving the directory structure.
        No zipping is performed - files are stored as-is on the server.

        Parameters
        ----------
        address : str
            DeepchemAddress
        dir_path : str
            Path to directory on disk
        card : dict, optional
            Metadata card
        kind : str
            Type of object

        Returns
        -------
        str
            The assigned address
        """
        profile, project, key = self._parse_address(address)

        # Collect all files with their relative paths
        files_to_upload = []
        for root, dirs, files in os.walk(dir_path):
            for file in files:
                file_path = os.path.join(root, file)
                relative_path = os.path.relpath(file_path, dir_path)
                with open(file_path, "rb") as f:
                    files_to_upload.append(("files", (relative_path, f.read())))

        # Build form data
        form_data: Dict[str, Any] = {"kind": kind}
        if card:
            form_data["card"] = json.dumps(card)

        response = self._get_client().post(
            f"{self.base_url}/directory/{profile}/{project}/{key}",
            files=files_to_upload,
            data=form_data,
        )
        response.raise_for_status()
        return response.json()["address"]

    def download_directory(
        self,
        address: str,
        dest_dir: str,
    ) -> str:
        """Download a directory by fetching files individually.

        This is used for internal operations (model loading, dataset loading).
        Files are downloaded individually and reassembled locally.

        Parameters
        ----------
        address : str
            DeepchemAddress
        dest_dir : str
            Destination directory on disk

        Returns
        -------
        str
            Path to destination directory
        """
        profile, project, key = self._parse_address(address)

        # Get list of files in the directory
        response = self._get_client().get(f"{self.base_url}/directory-contents/{profile}/{project}/{key}")
        response.raise_for_status()
        files = response.json()["files"]

        # Create destination directory
        os.makedirs(dest_dir, exist_ok=True)

        # Download each file
        for relative_path in files:
            file_key = f"{key}/{relative_path}"
            data = self.get_data(self._make_address(profile, project, file_key))

            dest_path = os.path.join(dest_dir, relative_path)
            os.makedirs(os.path.dirname(dest_path) or ".", exist_ok=True)
            with open(dest_path, "wb") as f:
                f.write(data)

        return dest_dir

    def download_directory_as_zip(
        self,
        address: str,
        dest_path: str,
    ) -> str:
        """Download a directory as a ZIP archive.

        This is for user-facing downloads where a single ZIP file is preferred.
        The ZIP is created on-demand by the server.

        Parameters
        ----------
        address : str
            DeepchemAddress
        dest_path : str
            Destination path for the ZIP file

        Returns
        -------
        str
            Path to downloaded ZIP file
        """
        profile, project, key = self._parse_address(address)

        response = self._get_client().get(
            f"{self.base_url}/download/{profile}/{project}/{key}",
            params={"format": "zip"},
        )
        response.raise_for_status()

        # Ensure parent directory exists
        os.makedirs(os.path.dirname(dest_path) or ".", exist_ok=True)

        with open(dest_path, "wb") as f:
            f.write(response.content)

        return dest_path

    def list_directory_contents(
        self,
        address: str,
    ) -> List[str]:
        """List all files in a directory.

        Parameters
        ----------
        address : str
            DeepchemAddress

        Returns
        -------
        list of str
            List of relative file paths
        """
        profile, project, key = self._parse_address(address)

        response = self._get_client().get(f"{self.base_url}/directory-contents/{profile}/{project}/{key}")
        response.raise_for_status()
        return response.json()["files"]

    def upload_from_profile_project(
        self,
        profile: str,
        project: str,
        key: str,
        data: bytes,
        card: Optional[Dict[str, Any]] = None,
        kind: str = "data",
    ) -> str:
        """Upload data using profile/project/key instead of address."""
        address = self._make_address(profile, project, key)
        return self.upload_data(address, data, card, kind)

    def get_from_profile_project(
        self,
        profile: str,
        project: str,
        key: str,
        fetch_sample: bool = False,
    ) -> bytes:
        """Get data using profile/project/key instead of address."""
        address = self._make_address(profile, project, key)
        return self.get_data(address, fetch_sample)
