import logging
import os
import shutil
import tempfile
from typing import Dict

from fastapi import APIRouter, BackgroundTasks, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse

from deepchem_server.core.cards import DataCard
from deepchem_server.utils import _download_data, _upload_data


logger = logging.getLogger("backend_logs")
logger.setLevel(logging.INFO)

router = APIRouter(prefix="/data", tags=["data"])

media_type_map = {
    ".csv": "text/csv",
    ".json": "application/json",
    ".txt": "text/plain",
    ".pdf": "application/pdf",
    ".png": "image/png",
    ".jpg": "image/jpeg",
    ".jpeg": "image/jpeg",
    ".zip": "application/zip",
    ".pdb": "chemical/x-pdb",
    ".sdf": "chemical/x-mdl-sdfile",
    ".xml": "application/xml",
}


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


@router.post("/download/data")
async def download_data(
        background_tasks: BackgroundTasks,
        profile_name: str = Form(...),
        project_name: str = Form(...),
        data_address: str = Form(...),
        backend: str = Form("local"),
) -> FileResponse:
    """
    Download data from datastore

    Parameters
    ----------
    background_tasks: BackgroundTasks
        FastAPI background tasks for cleanup
    profile_name: str
        Name of the Profile where the data is stored
    project_name: str
        Name of the Project where the data is stored
    data_address: str
        DeepChem address of the data object to download
    backend: str
        Backend to be used (Default: local)

    Returns
    -------
    FileResponse
        The file or zipped directory for download

    Raises
    ------
    HTTPException
        404 if the object doesn't exist
        400 if the address is invalid
        500 for other errors
    """
    try:
        file_path, is_directory, object_name = _download_data(
            profile_name=profile_name,
            project_name=project_name,
            address=data_address,
            backend=backend,
        )
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to download data: {str(e)}")

    if is_directory:
        temp_dir = tempfile.gettempdir()
        zip_base = tempfile.mktemp(dir=temp_dir, prefix=f"{object_name}_")

        try:
            zip_path = shutil.make_archive(zip_base, "zip", file_path)

            background_tasks.add_task(os.remove, zip_path)

            return FileResponse(path=zip_path, media_type="application/zip", filename=f"{object_name}.zip")
        except Exception as e:
            if os.path.exists(zip_base + ".zip"):
                os.remove(zip_base + ".zip")
            raise HTTPException(status_code=500, detail=f"Failed to create zip file: {str(e)}")
    else:
        _, ext = os.path.splitext(object_name)
        media_type = media_type_map.get(ext.lower(), "application/octet-stream")

        return FileResponse(path=file_path, media_type=media_type, filename=object_name)


@router.post("/download/model")
async def download_model(
        background_tasks: BackgroundTasks,
        profile_name: str = Form(...),
        project_name: str = Form(...),
        model_address: str = Form(...),
        backend: str = Form("local"),
) -> FileResponse:
    """
    Download model from datastore

    Parameters
    ----------
    background_tasks: BackgroundTasks
        FastAPI background tasks for cleanup
    profile_name: str
        Name of the Profile where the model is stored
    project_name: str
        Name of the Project where the model is stored
    model_address: str
        DeepChem address of the model object to download
    backend: str
        Backend to be used (Default: local)

    Returns
    -------
    FileResponse
        The zipped model directory for download

    Raises
    ------
    HTTPException
        404 if the object doesn't exist
        400 if the address is invalid
        500 for other errors
    """
    try:
        file_path, is_directory, object_name = _download_data(
            profile_name=profile_name,
            project_name=project_name,
            address=model_address,
            backend=backend,
        )
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to download model: {str(e)}")

    # Models are typically directories, but handle both cases
    if is_directory:
        # Create a temporary zip file
        temp_dir = tempfile.gettempdir()
        # Use a unique temp name without extension, make_archive adds .zip
        zip_base = tempfile.mktemp(dir=temp_dir, prefix=f"{object_name}_")

        try:
            # Create zip archive (make_archive adds .zip extension automatically)
            zip_path = shutil.make_archive(zip_base, "zip", file_path)

            # Schedule cleanup of the zip file after response is sent
            background_tasks.add_task(os.remove, zip_path)

            return FileResponse(path=zip_path, media_type="application/zip", filename=f"{object_name}.zip")
        except Exception as e:
            # Clean up zip file if it was created
            if os.path.exists(zip_base + ".zip"):
                os.remove(zip_base + ".zip")
            raise HTTPException(status_code=500, detail=f"Failed to create zip file: {str(e)}")
    else:
        # For single file models (unlikely but handle it)
        return FileResponse(path=file_path, media_type="application/octet-stream", filename=object_name)
