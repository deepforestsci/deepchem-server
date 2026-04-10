"""
Data client module for interacting with DeepChem Server data endpoints.

Contains the Data class for all data management operations.
"""
import io
import json
import tempfile
import zipfile
from pathlib import Path
from typing import Any, Dict, Optional, Union

from requests_toolbelt import MultipartEncoder

from deepchem_server.core.common.address import DeepchemAddress

from .base import BaseClient
from .settings import Settings


class Data(BaseClient):
    """
    Client for interacting with DeepChem Server data endpoints.

    This class provides methods for data management operations.
    """

    def __init__(self, settings: Optional[Settings] = None, base_url: Optional[str] = None):
        """
        Initialize Data client.

        Parameters
        ----------
            settings: Settings
                Settings instance for configuration
            base_url: str
                Base URL for the API (overrides settings if provided)
        """
        super().__init__(settings, base_url)

    def upload_data(
        self,
        file_path: Union[str, Path],
        filename: Optional[str] = None,
        description: Optional[str] = None,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
        backend: str = "local",
    ) -> Dict[str, Any]:
        """
        Upload data to datastore.

        Args:
            file_path: Path to the file to upload
            filename: File name to save the uploaded file (uses original name if not provided)
            description: Description of the file
            profile_name: Profile name (uses settings if not provided)
            project_name: Project name (uses settings if not provided)
            backend: Backend to be used to run the job (Default: local)

        Returns:
            Response containing the dataset address

        Raises:
            ValueError: If required settings are missing or file doesn't exist
            requests.exceptions.RequestException: If API request fails
        """
        file_path = Path(file_path)

        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        # Get profile and project names (validates configuration)
        profile, project = self._get_profile_and_project(profile_name, project_name)

        if filename is None:
            filename = file_path.name

        # Prepare the multipart form data using MultipartEncoder
        fields: dict[str, Any] = {
            "profile_name": profile,
            "project_name": project,
            "filename": filename,
            "backend": backend,
            "file": (filename, open(file_path, "rb")),
        }

        if description is not None:
            fields["description"] = description

        multipart_data = MultipartEncoder(fields=fields)

        try:
            response = self._post(
                "/data/uploaddata",
                data=multipart_data,
                headers={"Content-Type": multipart_data.content_type},
            )
            fields["file"][1].close()
            return self._validate_response(response)
        except Exception as e:
            fields["file"][1].close()
            raise ValueError(f"Failed to upload data: {str(e)}") from e

    def list_data(self, profile_name: Optional[str] = None, project_name: Optional[str] = None) -> Dict[str, Any]:
        """
        List data in datastore.

        Parameters
        ----------
            profile_name: str
                Profile name (uses settings if not provided)
            project_name: str
                Project name (uses settings if not provided)

        Returns
        -------
            Dict[str, Any]
                List of data in datastore

        Raises
        ------
            requests.exceptions.HTTPError: If response indicates an error

        Examples
        --------
        >>> from pyds.data import Data
        >>> client = Data()
        >>> response = client.list_data()
        >>> print(response)
        ["data1.csv", "data2.csv"]
        """
        profile, project = self._get_profile_and_project(profile_name, project_name)
        response = self._get(f"/data?profile_name={profile}&project_name={project}")
        data = self._validate_response(response)
        return data["data_files"]

    def download_data(
        self,
        address: str,
        destination_path: str,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> str:
        """
        Download data from datastore.

        Parameters
        ----------
            address: str
                Address of the data to download
            destination_path: str
                Path to save the downloaded file
            profile_name: str
                Profile name (uses settings if not provided)
            project_name: str
                Project name (uses settings if not provided)

        Returns
        -------
            str
                Absolute path to the downloaded file

        Raises
        ------
            requests.exceptions.HTTPError: If response indicates an error

        Examples
        --------
        >>> from pyds.data import Data
        >>> client = Data()
        >>> response = client.download_data("data1.csv", "/path/to/download/data.csv")
        >>> print(response)
        "/path/to/download/data.csv"

        """
        profile, project = self._get_profile_and_project(profile_name, project_name)
        response = self._get(f"/data/{address}?profile_name={profile}&project_name={project}")
        return self._validate_file_response(response, destination_path)

    def get(
        self,
        address: str,
        destination_path: Optional[str] = None,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Any:
        """
        Download and parse data from the datastore.

        The return type is inferred from the file extension:

        Data -> Python object
        .csv -> pandas.DataFrame
        .json -> dict / list
        .png -> PIL.Image.Image
        .txt, .pdb, .fasta -> str
        .fastq, .sdf, .xml -> str
        .pdbqt, .smi etc. -> str
        .zip / dataset directory -> dc.data.DiskDataset
        anything else -> bytes

        If destination_path is supplied the raw file is also saved
        there.  When omitted a named temporary file is used for formats
        that require disk access (e.g. DiskDataset).

        Parameters
        ----------
        address : str
            Address of the data to download.
        destination_path : str, optional
            Path to save the raw file. Auto-generated when not provided.
        profile_name : str, optional
            Profile name (uses settings if not provided).
        project_name : str, optional
            Project name (uses settings if not provided).

        Returns
        -------
        Any
            Parsed data object.
        """
        profile, project = self._get_profile_and_project(profile_name, project_name)
        key = DeepchemAddress.get_key(address)
        response = self._get(f"/data/{key}?profile_name={profile}&project_name={project}")
        self._handle_http_error(response)

        if destination_path is not None:
            with open(destination_path, "wb") as f:
                f.write(response.content)

        parsed_data = self._parse_content(response.content, key, response.headers.get("content-type", ""))
        if parsed_data is None:
            return f"Unkown data type - Unable to get data {address}"
        return parsed_data

    _TEXT_EXTENSIONS = frozenset({
        "pdb",
        "fasta",
        "fastq",
        "sdf",
        "txt",
        "xml",
        "pdbqt",
        "smi",
        "smiles",
        "cxsmiles",
        "py",
        "log",
    })

    def _parse_content(self, content: bytes, address: str, content_type: str = "") -> Any:
        """Parse raw bytes into a Python object based on the file extension."""
        ext = address.rsplit(".", 1)[-1].lower() if "." in address else ""

        if ext == "csv":
            import pandas as pd

            return pd.read_csv(io.BytesIO(content))

        if ext == "json":
            return json.loads(content.decode("utf-8"))

        if ext == "png":
            try:
                from PIL import Image

                return Image.open(io.BytesIO(content))
            except ImportError:
                return content

        if ext in self._TEXT_EXTENSIONS:
            return content.decode("utf-8", errors="replace")

        if ext == "zip" or content_type == "application/zip":
            extract_dir = tempfile.mkdtemp(prefix="dc_dataset_")
            with zipfile.ZipFile(io.BytesIO(content)) as zf:
                zf.extractall(extract_dir)
            try:
                import deepchem as dc
                import os

                entries = os.listdir(extract_dir)
                if len(entries) == 1:
                    potential_dir = os.path.join(extract_dir, entries[0])
                    if os.path.isdir(potential_dir):
                        extract_dir = potential_dir

                return dc.data.DiskDataset(data_dir=extract_dir)
            except ImportError:
                return {
                    "path": extract_dir,
                    "note": "Install deepchem to load as dc.data.DiskDataset.",
                }

        return None
