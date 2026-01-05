"""
Data Router - Endpoints for data upload and download operations.
"""
import logging
import os
from typing import Dict

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import StreamingResponse

from deepchem_server.api.utils import upload_data as _upload_data
from deepchem_server.core.common.cards import DataCard

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

    file_type = filename.split('.')[-1]  # getting extension
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


@router.get("/download")
async def download_data(
        profile_name: str,
        project_name: str,
        address: str,
        format: str = "auto",
):
    """
    Download data from datastore

    Parameters
    ----------
    profile_name: str
        Name of the Profile
    project_name: str
        Name of the Project
    address: str
        DeepChem address of the object to download
    format: str
        Download format: 'zip' for directories, 'auto' for automatic detection

    Returns
    -------
    dict
        File content or redirect to datastore download endpoint

    Notes
    -----
    For directories, format='zip' will return a ZIP archive.
    For files, the file content is returned directly.
    """
    import httpx

    from deepchem_server.core import config

    datastore = config.get_datastore()
    if datastore is None:
        raise HTTPException(status_code=500, detail="Datastore not configured")

    # Build the download URL for the datastore service
    datastore_url = os.getenv("DATASTORE_URL", "http://localhost:8081")
    api_key = os.getenv("DATASTORE_API_KEY", "dev-api-key")

    # Parse the address to get profile/project/key
    if address.startswith("deepchem://"):
        address_path = address[len("deepchem://"):]
    else:
        address_path = f"{profile_name}/{project_name}/{address}"

    download_url = f"{datastore_url}/api/v1/download/{address_path}"

    async def stream_response():
        async with httpx.AsyncClient() as client:
            async with client.stream(
                "GET",
                download_url,
                params={"format": format},
                headers={"X-API-Key": api_key},
            ) as response:
                async for chunk in response.aiter_bytes():
                    yield chunk

    # Get headers from datastore for content-disposition
    async with httpx.AsyncClient() as client:
        head_response = await client.head(
            download_url,
            params={"format": format},
            headers={"X-API-Key": api_key},
        )

    content_type = head_response.headers.get("content-type", "application/octet-stream")
    content_disposition = head_response.headers.get("content-disposition", "")

    return StreamingResponse(
        stream_response(),
        media_type=content_type,
        headers={"Content-Disposition": content_disposition} if content_disposition else {},
    )
