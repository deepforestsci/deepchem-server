import logging
from typing import Dict, Optional
import mimetypes
import os

from fastapi import APIRouter, File, Form, UploadFile, HTTPException
from fastapi.responses import FileResponse
from starlette.background import BackgroundTask

from deepchem_server.utils import _init_datastore
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.common.address import DeepchemAddress
from deepchem_server.utils import _upload_data, _download_file


logger = logging.getLogger("backend_logs")
logger.setLevel(logging.INFO)

router = APIRouter(prefix="/data", tags=["data"])


@router.post("/uploaddata")
async def upload_data(
        file: UploadFile = File(...),
        profile_name: str = Form(...),
        project_name: str = Form(...),
        filename: str = Form(None),
        description: str = Form(None),
        backend="local",
) -> Dict:
    """
    Upload data to datastore

    Parameters
    ----------
    file: UploadFile
        A file uploaded in a request
    profile_name: str
        Name of the Profile where the job is run
    project_name: str
        Name of the Project where the job is run
    filename: str
        File name to save the uploaded file
    description: Union[Dict, str]
        Description of the file
    backend: str
        Backend to be used to run the job (Default: local)
    """
    contents = await file.read()

    if filename is None:
        filename = file.filename

    file_type = filename.split('.')[-1]
    if file_type in ['csv', 'parquet']:
        data_type = 'pandas.DataFrame'
    elif file_type in [
            "pdb",
            "sdf",
            "fasta",
            "fastq",
            "sdf",
            "txt",
            "xml",
            "pdbqt",
            "smi",
            "smiles",
            "cxsmiles",
            "json",
    ]:
        data_type = 'text/plain'
    elif file_type in ['dcd', 'bz2', 'zip', 'onnx', 'hdf5']:
        data_type = 'binary'
    elif file_type in ['png']:
        data_type = 'png'
    else:
        data_type = ''

    card: DataCard = DataCard(address='', file_type=file_type, data_type=data_type, description=description)

    address: str = _upload_data(profile_name, project_name, filename, contents, card, backend=backend)  # type: ignore
    return {"dataset_address": address}


@router.get("")
async def list_files(profile_name: str, project_name: str, backend="local") -> Dict:
    """
    List data files in datastore for a given profile and project

    Parameters
    ----------
    profile_name: str
        Name of the Profile where the job is run
    project_name: str
        Name of the Project where the job is run
    backend: str
        Backend to be used to run the job (Default: local)
    """
    datastore = _init_datastore(profile_name=profile_name, project_name=project_name, backend=backend)
    data_files = datastore.list_data().split("\n")
    return {"data_files": data_files}


@router.get("/{file_name:path}")
async def download_file(profile_name: Optional[str],
                        project_name: Optional[str],
                        file_name: str,
                        backend="local") -> FileResponse:
    """
    Download a file from the datastore for a given profile and project

    Parameters
    ----------
    profile_name: str
        Name of the Profile where the job is run
    project_name: str
        Name of the Project where the job is run
    file_name: str
        Name of the file to download
    backend: str
        Backend to be used to run the job (Default: local)

    Returns
    -------
    FileResponse:
        The file response object

    Raises
    ------
    HTTPException:
        If the file is not found in the datastore
    """
    if file_name.startswith(DeepchemAddress.address_prefix):
        parsed_address = DeepchemAddress(file_name).parse_address(file_name)
        profile_name = parsed_address["profile"]
        project_name = parsed_address["project"]
        file_name = parsed_address["key"]

    if profile_name is None or project_name is None:
        raise HTTPException(status_code=400, detail="Profile and project names are required")

    address = f"deepchem://{profile_name}/{project_name}/{file_name}"
    print(address)
    key = DeepchemAddress(address).key
    file_path = _download_file(profile_name=profile_name,
                               project_name=project_name,
                               file_name=file_name,
                               backend=backend)
    if file_path is None:
        raise HTTPException(status_code=404, detail=f"The address {file_name} is not found in the datastore")
    mime_type, _ = mimetypes.guess_type(file_path)
    if mime_type is None:
        mime_type = "application/octet-stream"

    def clean_up():
        if backend == "local":
            return

        if os.path.exists(file_path):
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                import shutil

                shutil.rmtree(file_path)
        if mime_type == "application/zip":
            zip_dir_path = os.path.join(os.path.dirname(file_path), key)
            if os.path.exists(zip_dir_path):
                import shutil

                shutil.rmtree(zip_dir_path)

    return FileResponse(file_path, filename=key, media_type=mime_type, background=BackgroundTask(clean_up))
