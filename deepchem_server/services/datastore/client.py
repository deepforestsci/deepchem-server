"""
Python client for DeepchemDatastore service.

Provides a convenient interface for other services to communicate
with the datastore service.
"""

import json
import logging
from typing import Any, Dict, List, Optional

import httpx


logger = logging.getLogger(__name__)


class DatastoreClient:
    """Client for communicating with DeepchemDatastore service.

    Parameters
    ----------
    host : str
        Hostname of the datastore service
    port : int
        Port of the datastore service
    api_key : str
        API key for authentication
    timeout : float
        Request timeout in seconds
    """

    def __init__(
        self,
        host: str = "localhost",
        port: int = 8081,
        api_key: str = "dev-api-key-change-in-production",
        timeout: float = 30.0,
    ) -> None:
        self.base_url = f"http://{host}:{port}/api/v1"
        self.api_key = api_key
        self.timeout = timeout
        self._client: Optional[httpx.Client] = None

    def _get_client(self) -> httpx.Client:
        """Get or create the HTTP client."""
        if self._client is None:
            self._client = httpx.Client(
                timeout=self.timeout,
                headers={"X-API-Key": self.api_key},
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
        form_data = {"kind": kind}
        if card:
            form_data["card"] = json.dumps(card)

        response = self._get_client().post(
            f"{self.base_url}/data/{profile}/{project}/{key}",
            files=files,
            data=form_data,
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
        """Delete an object from the datastore.

        Parameters
        ----------
        address : str
            DeepchemAddress
        kind : str
            Type of object

        Returns
        -------
        bool
            True if successful
        """
        profile, project, key = self._parse_address(address)

        response = self._get_client().delete(
            f"{self.base_url}/data/{profile}/{project}/{key}",
            params={"kind": kind},
        )
        response.raise_for_status()
        return response.json()["status"] == "success"

    def list_data(self, profile: str, project: str) -> List[str]:
        """List all keys in a profile/project.

        Parameters
        ----------
        profile : str
            Profile name
        project : str
            Project name

        Returns
        -------
        list of str
            List of keys
        """
        response = self._get_client().get(f"{self.base_url}/list/{profile}/{project}",)
        response.raise_for_status()
        return response.json()["keys"]

    def get_card(
        self,
        address: str,
        kind: str = "data",
    ) -> Optional[Dict[str, Any]]:
        """Get the metadata card for an object.

        Parameters
        ----------
        address : str
            DeepchemAddress
        kind : str
            Type of object

        Returns
        -------
        dict or None
            Card metadata
        """
        profile, project, key = self._parse_address(address)

        response = self._get_client().get(
            f"{self.base_url}/card/{profile}/{project}/{key}",
            params={"kind": kind},
        )
        response.raise_for_status()
        return response.json()["card"]

    def exists(self, address: str) -> bool:
        """Check if an object exists.

        Parameters
        ----------
        address : str
            DeepchemAddress

        Returns
        -------
        bool
            True if exists
        """
        profile, project, key = self._parse_address(address)

        response = self._get_client().get(f"{self.base_url}/exists/{profile}/{project}/{key}",)
        response.raise_for_status()
        return response.json()["exists"]

    def get_size(self, address: str) -> int:
        """Get the size of an object.

        Parameters
        ----------
        address : str
            DeepchemAddress

        Returns
        -------
        int
            Size in bytes
        """
        profile, project, key = self._parse_address(address)

        response = self._get_client().get(f"{self.base_url}/size/{profile}/{project}/{key}",)
        response.raise_for_status()
        return response.json()["size"]

    def healthcheck(self) -> bool:
        """Check if the service is healthy.

        Returns
        -------
        bool
            True if healthy
        """
        try:
            # Healthcheck doesn't require auth
            response = httpx.get(
                f"{self.base_url}/healthcheck",
                timeout=5.0,
            )
            return response.status_code == 200
        except Exception:
            return False
