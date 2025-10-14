import logging
from typing import Dict

from fastapi import APIRouter, File, Form, UploadFile

from deepchem_server.core.cards import DataCard
from deepchem_server.utils import _upload_data


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
