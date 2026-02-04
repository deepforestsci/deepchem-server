"""
DeepchemDatastore - High-level datastore interface for deepchem operations.

This class provides the interface that all deepchem primitives expect.
It communicates with the datastore service via HTTP and handles:
- File upload/download
- Model loading (downloads, extracts, instantiates)
- Dataset loading (downloads, extracts, creates DiskDataset)
- Checkpoint scanning for multicore featurization
"""
import io
import json
import logging
import os
import shutil
import tempfile
from typing import Any, BinaryIO, List, Optional, Union

from httpx import HTTPStatusError
import pandas as pd

# DeepChem imports - lazy loaded to avoid import errors if not installed
# but this package does require deepchem as a dependency for full functionality
try:
    import deepchem as dc
    HAS_DEEPCHEM = True
except ImportError:
    HAS_DEEPCHEM = False
    dc = None

# PIL imports for image handling
try:
    from PIL.PngImagePlugin import PngImageFile
    HAS_PIL = True
except ImportError:
    HAS_PIL = False
    PngImageFile = None  # type: ignore

from deepchem_server.core.common.address import DeepchemAddress
from deepchem_server.core.common.cards import Card, DataCard, ModelCard
from deepchem_server.services.datastore.client.http import DatastoreClient


logger = logging.getLogger(__name__)


def _get_csv_or_dataframe_shape(filename: Optional[str] = None, dataframe: Optional[pd.DataFrame] = None) -> tuple:
    """Get shape of CSV file or DataFrame."""
    if dataframe is not None:
        return (len(dataframe), len(dataframe.columns))
    elif filename is not None:
        df = pd.read_csv(filename)
        return (len(df), len(df.columns))
    return (0, 0)


class DeepchemDatastore:
    """High-level datastore interface for deepchem operations.
    
    This class provides the interface that all primitives (feat, train, evaluate)
    expect. It wraps the HTTP DatastoreClient and handles complex operations
    like model/dataset loading.
    
    Parameters
    ----------
    client : DatastoreClient
        The HTTP client for communicating with DatastoreService
    profile : str
        Profile name (used for address construction)
    project : str
        Project name (used for address construction)
    temp_dir : str, optional
        Directory for temporary file operations (downloaded models/datasets)
    
    Attributes
    ----------
    storage_loc : str
        Local temp directory for downloaded files (for compatibility)
    address_prefix : str
        Prefix for constructing addresses (e.g., "profile/project/")
    sample_rows : int
        Number of rows to return when fetch_sample=True (default: 100)
    
    Examples
    --------
    >>> from deepchem_server.services.datastore.client import DatastoreClient, DeepchemDatastore
    >>> client = DatastoreClient(url="http://localhost:8081")
    >>> datastore = DeepchemDatastore(client, "demo_user", "demo_project")
    >>>
    >>> # Upload data
    >>> address = datastore.upload_data("data.csv", "/path/to/data.csv", card)
    >>>
    >>> # Get data
    >>> df = datastore.get("deepchem://demo_user/demo_project/data.csv")
    """

    def __init__(
        self,
        client: DatastoreClient,
        profile_name: str,
        project_name: str,
        basedir: Optional[str] = None,
    ) -> None:
        self.client = client
        self.profile = profile_name
        self.project = project_name
        self.address_prefix = f"{profile_name}/{project_name}/"
        self.storage_loc = basedir or tempfile.mkdtemp(prefix="datastore_")
        os.makedirs(self.storage_loc, exist_ok=True)
        self.sample_rows = 100

    def _make_address(self, key: str) -> str:
        """Create a deepchem address from a key."""
        if key.startswith("deepchem://"):
            return key
        return f"deepchem://{self.profile}/{self.project}/{key}"

    def upload_data(
        self,
        datastore_filename: str,
        filename: str,
        card: Union[ModelCard, DataCard],
        kind: Optional[str] = "data",
    ) -> str:
        """Upload data to the datastore.
        
        Parameters
        ----------
        datastore_filename : str
            The name within the datastore
        filename : str
            Path to file/directory on disk
        card : ModelCard or DataCard
            Metadata card
        kind : str, optional
            Type of object
            
        Returns
        -------
        str
            The assigned address
        """
        address = self._make_address(datastore_filename)
        card.address = address

        # Set shape for CSV files
        if isinstance(card, DataCard) and datastore_filename.endswith('.csv'):
            card.shape = _get_csv_or_dataframe_shape(filename=filename)

        card_dict = json.loads(card.to_json()) if card else None

        if os.path.isdir(filename):
            # Directory: ZIP and upload
            return self.client.upload_directory(address, filename, card_dict, kind or "data")
        else:
            # File: upload directly
            return self.client.upload_file(address, filename, card_dict, kind or "data")

    def upload_data_from_memory(
        self,
        data: Any,
        datastore_filename: str,
        card: Union[DataCard, ModelCard, None],
        kind: str = "data",
    ) -> str:
        """Upload in-memory data to the datastore.
        
        Supports the following data types:
        - Card: Serialized card data
        - pd.DataFrame: CSV format
        - dc.data.NumpyDataset: Converted to DiskDataset and uploaded as ZIP
        - dc.data.DiskDataset: Uploaded as ZIP
        - dc.models.Model: Model directory uploaded as ZIP
        - str: Text file
        - bytes: Binary file
        - PIL.PngImagePlugin.PngImageFile: PNG image
        
        Parameters
        ----------
        data : Any
            Data to upload (DataFrame, DiskDataset, Model, bytes, str, etc.)
        datastore_filename : str
            The name within the datastore
        card : DataCard, ModelCard, or None
            Metadata card
        kind : str
            Type of object ('data' or 'model')
            
        Returns
        -------
        str
            The assigned address
            
        Raises
        ------
        ValueError
            If unsupported data type is provided
        """
        address = self._make_address(datastore_filename)

        if card is not None:
            card.address = address

        # Handle Card objects (just serialize)
        if isinstance(data, Card):
            data_bytes = bytes(data)
            card_dict = json.loads(card.to_json()) if card else None
            return self.client.upload_data(address, data_bytes, card_dict, kind)

        # Handle pandas DataFrame
        if isinstance(data, pd.DataFrame):
            buffer = io.BytesIO()
            data.to_csv(buffer, index=False)
            data_bytes = buffer.getvalue()
            if isinstance(card, DataCard):
                card.shape = (len(data), len(data.columns))
            card_dict = json.loads(card.to_json()) if card else None
            return self.client.upload_data(address, data_bytes, card_dict, kind)

        # Handle DeepChem NumpyDataset - convert to DiskDataset first
        if HAS_DEEPCHEM and isinstance(data, dc.data.NumpyDataset):
            if isinstance(card, DataCard):
                card.shape = data.get_shape()
            # Create temporary DiskDataset
            temp_dir = tempfile.mkdtemp(prefix="numpy_to_disk_")
            dest_path = os.path.join(temp_dir, datastore_filename)
            dc.data.DiskDataset.from_numpy(data.X, data.y, data.w, data.ids, data_dir=dest_path)
            card_dict = json.loads(card.to_json()) if card else None
            try:
                return self.client.upload_directory(address, dest_path, card_dict, kind)
            finally:
                shutil.rmtree(temp_dir, ignore_errors=True)

        # Handle DeepChem DiskDataset
        if HAS_DEEPCHEM and isinstance(data, dc.data.DiskDataset):
            if isinstance(card, DataCard):
                card.shape = data.get_shape()
            card_dict = json.loads(card.to_json()) if card else None
            return self.client.upload_directory(address, data.data_dir, card_dict, kind)

        # Handle DeepChem Model
        if HAS_DEEPCHEM and hasattr(data, 'model_dir') and isinstance(getattr(data, 'model_dir', None), str):
            # This is likely a DeepChem model
            card_dict = json.loads(card.to_json()) if card else None
            return self.client.upload_directory(address, data.model_dir, card_dict, kind)

        # Handle PIL PngImageFile
        if HAS_PIL and PngImageFile is not None and isinstance(data, PngImageFile):
            buffer = io.BytesIO()
            data.save(buffer, format='PNG')
            data_bytes = buffer.getvalue()
            card_dict = json.loads(card.to_json()) if card else None
            return self.client.upload_data(self._make_address(address), data_bytes, card_dict, kind)

        # Handle string (text file)
        if isinstance(data, str):
            data_bytes = data.encode("utf-8")
            card_dict = json.loads(card.to_json()) if card else None
            return self.client.upload_data(self._make_address(address), data_bytes, card_dict, kind)

        # Handle bytes (binary file)
        if isinstance(data, bytes):
            card_dict = json.loads(card.to_json()) if card else None
            return self.client.upload_data(self._make_address(address), data, card_dict, kind)

        # Unsupported type
        raise ValueError(f"Unsupported data type: {type(data)}. "
                         f"Supported types: Card, pd.DataFrame, dc.data.NumpyDataset, "
                         f"dc.data.DiskDataset, dc.models.Model, str, bytes, PngImageFile")

    def upload_model(
        self,
        modelname: str,
        model: Any,
        card: ModelCard,
    ) -> str:
        """Upload a model to the datastore.
        
        Parameters
        ----------
        modelname : str
            Name of the model
        model : dc.models.Model
            DeepChem model to upload
        card : ModelCard
            Model metadata
            
        Returns
        -------
        str
            The assigned address
        """
        address = self._make_address(modelname)
        card.address = address
        card_dict = json.loads(card.to_json())

        # Models are directories - ZIP and upload
        return self.client.upload_directory(self._make_address(address), model.model_dir, card_dict, "model")

    def upload_model_from_memory(
        self,
        model_name: str,
        model_files: List[BinaryIO],
        model_filenames: List[str],
        card: ModelCard,
    ) -> str:
        """Upload model files from memory to the datastore.
        
        Parameters
        ----------
        model_name : str
            Name of the model
        model_files : list of file-like objects
            List of open file handles for model files
        model_filenames : list of str
            List of filenames for each file
        card : ModelCard
            Model metadata
            
        Returns
        -------
        str
            The assigned address
        """
        address = self._make_address(model_name)
        card.address = address

        # Create a temporary directory to store files
        temp_dir = tempfile.mkdtemp(prefix="model_upload_")
        model_dir = os.path.join(temp_dir, model_name)
        os.makedirs(model_dir, exist_ok=True)

        try:
            # Write all files to temp directory
            for file_handle, filename in zip(model_files, model_filenames):
                file_path = os.path.join(model_dir, filename)
                os.makedirs(os.path.dirname(file_path), exist_ok=True)
                with open(file_path, 'wb') as f:
                    f.write(file_handle.read())

            card_dict = json.loads(card.to_json())
            return self.client.upload_directory(self._make_address(address), model_dir, card_dict, "model")
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

    def get(
        self,
        address: str,
        kind: Optional[str] = "data",
        fetch_sample: bool = False,
    ) -> Any:
        """Fetch something from datastore at address.
        
        Parameters
        ----------
        address : str
            DeepchemAddress
        kind : str, optional
            'data' or 'model'
        fetch_sample : bool
            Whether to get sample or full data
            
        Returns
        -------
        Any
            The requested object
        """
        if address.endswith(".cdc"):
            return self.get_card(address[:-4], kind="data")
        elif address.endswith(".cmc"):
            return self.get_card(address[:-4], kind="model")

        if kind == "data":
            return self.get_data(address, fetch_sample)
        elif kind == "model":
            return self.get_model(address)
        return None

    def get_data(self, address: str, fetch_sample: bool = False) -> Any:
        """Fetch data from datastore.
        
        Parameters
        ----------
        address : str
            DeepchemAddress
        fetch_sample : bool
            Whether to fetch only a sample
            
        Returns
        -------
        Any
            The data object (DataFrame, DiskDataset, etc.)
        """
        card = self.get_card(address, kind="data")
        if card and isinstance(card, DataCard):
            if card.file_type == "dir" or card.data_type == "dc.data.DiskDataset":
                if not HAS_DEEPCHEM:
                    raise ImportError("deepchem is required to load DiskDataset")
                key = DeepchemAddress.get_key(address)
                dest_dir = os.path.join(self.storage_loc, key)

                try:
                    self.client.download_directory(self._make_address(address), dest_dir)
                except HTTPStatusError as e:
                    if e.response.status_code == 404:
                        raise FileNotFoundError(f"Object not found: {address}")
                    else:
                        raise ValueError(f"Error downloading directory: {e}")
                except Exception as e:
                    raise ValueError(f"Error downloading directory: {e}")

                return dc.data.DiskDataset(data_dir=dest_dir)

            else:
                try:
                    data_bytes = self.client.get_data(self._make_address(address), fetch_sample)
                except HTTPStatusError as e:
                    if e.response.status_code == 404:
                        raise FileNotFoundError(f"Object not found: {address}")
                    else:
                        raise ValueError(f"Error downloading file: {e}")
                except Exception as e:
                    raise ValueError(f"Error downloading file: {e}")

                if card.file_type == "csv":
                    if fetch_sample:
                        return pd.read_csv(io.BytesIO(data_bytes), nrows=self.sample_rows)
                    return pd.read_csv(io.BytesIO(data_bytes))
                elif card.file_type == "json":
                    return json.loads(data_bytes.decode("utf-8"))
                elif card.file_type == "txt":
                    return data_bytes.decode("utf-8").splitlines(keepends=True)
                elif card.file_type == "xml":
                    return data_bytes.decode("utf-8").splitlines(keepends=True)
                elif card.file_type == "png":
                    if HAS_PIL:
                        from PIL import Image

                        return Image.open(io.BytesIO(data_bytes))
                    return data_bytes

                # Default: return raw bytes
                return data_bytes

        raise FileNotFoundError(f"Object not found: {address}")

    def get_model(self, address: str) -> Any:
        """Fetch and load a model from datastore.
        
        Parameters
        ----------
        address : str
            DeepchemAddress
            
        Returns
        -------
        Any
            The loaded DeepChem model
        """
        from deepchem_server.core.common import model_mappings

        card = self.get_card(address, kind="model")
        if not isinstance(card, ModelCard):
            raise ValueError(f"Expected ModelCard for address {address}, got {type(card)}")

        # Download model directory
        key = DeepchemAddress.get_key(address)
        dest_dir = os.path.join(self.storage_loc, key)
        self.client.download_directory(self._make_address(address), dest_dir)

        # Load model
        model = model_mappings.model_address_map[card.model_type](model_dir=dest_dir, **card.init_kwargs)
        try:
            model.restore()
        except AttributeError:
            model.reload()

        return model

    def get_card(
        self,
        address: str,
        kind: Optional[str] = "data",
    ) -> Optional[Union[DataCard, ModelCard]]:
        """Fetch card from datastore.
        
        Parameters
        ----------
        address : str
            DeepchemAddress
        kind : str, optional
            'data' or 'model'
            
        Returns
        -------
        DataCard, ModelCard, or None
            The card object
        """
        card_dict = self.client.get_card(self._make_address(address), kind or "data")

        if card_dict is None:
            return None

        if kind == "data":
            return DataCard.from_json(json.dumps(card_dict))
        elif kind == "model":
            return ModelCard.from_json(json.dumps(card_dict))
        return None

    def get_dir(self, address: str) -> str:
        """Download directory and return local path.
        
        Parameters
        ----------
        address : str
            DeepchemAddress
            
        Returns
        -------
        str
            Local path to downloaded directory
        """
        key = DeepchemAddress.get_key(address)
        dest_dir = os.path.join(self.storage_loc, key)
        self.client.download_directory(self._make_address(address), dest_dir)
        return dest_dir

    def exists(self, address: str) -> bool:
        """Check if an object exists."""
        return self.client.exists(self._make_address(address))

    def get_file_size(self, address: str) -> int:
        """Get size of an object."""
        return self.client.get_size(self._make_address(address))

    def list_data(self, include_card_files: bool = False) -> List[str]:
        """List all data in the datastore.

        Parameters
        ----------
        include_card_files : bool
            Whether to include card files (.cdc, .cmc)

        Returns
        -------
        List[str]
            List of data keys
        """
        return self.client.list_data(self.profile, self.project, include_card_files)

    def delete_object(self, address: str, kind: str = "data") -> bool:
        """Delete an object from datastore.
        
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
        return self.client.delete_object(self._make_address(address), kind)

    def download_object(
        self,
        address: str,
        filename: Optional[str] = None,
    ) -> None:
        """Download an object to disk.
        
        Parameters
        ----------
        address : str
            DeepchemAddress
        filename : str, optional
            Destination path
        """
        if not filename:
            raise ValueError("filename must be specified")

        card = self.get_card(self._make_address(address))
        if card and hasattr(card, "data_type") and card.data_type == "dc.data.DiskDataset":
            self.client.download_directory(self._make_address(address), filename)
        else:
            self.client.download_file(self._make_address(address), filename)

    def _get_datastore_objects(self, directory: str) -> List[str]:
        """Walk datastore and collect all objects as paths.
        
        This method is used by feat.py for checkpoint scanning during
        multicore featurization restart.
        
        Parameters
        ----------
        directory : str
            The directory to scan (typically self.storage_loc)
            
        Returns
        -------
        list of str
            List of full paths including card files (.cdc, .cmc)
        """
        return self.client.list_all_objects(self.profile, self.project, prefix=directory)

    def add_dir(self, dir_name: str) -> None:
        """Create a directory marker in the datastore.

        Parameters
        ----------
        dir_name : str
            The name of the directory to create

        Returns
        -------
        str
            The assigned address
        """
        try:
            self.client.create_directory(self._make_address(dir_name))
        except Exception as e:
            raise ValueError(f"Error creating directory: {e}")

    def move_object(self, address: str, destination: str) -> str:
        """Move an object to a new location.

        Parameters
        ----------
        address : str
            DeepchemAddress
        destination : str
            New location address
            If only the key is provided, then the client will assume
            the profile and project are the same as the original address,
            so this method can be assumed as mv command in the shell.

        Returns
        -------
        str
            The assigned address
        """
        return self.client.move_object(self._make_address(address), self._make_address(destination))
