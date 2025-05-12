import logging
from typing import Union, Dict
from fastapi import APIRouter, UploadFile, File, Form
from deepchem_server.core.cards import DataCard
from deepchem_server.utils import _upload_data


logger = logging.getLogger("backend_logs")
logger.setLevel(logging.INFO)

router = APIRouter(prefix="/data",
                   tags=["data"])


@router.post("/uploaddata")
async def upload_data(file: UploadFile = File(...),
                      profile_name: str = Form(None),
                      project_name: str = Form(None),
                      filename: str = Form(None),
                      description: Union[Dict, str] = Form(None),
                      backend='local'):
    """
    Upload data to datastore
    """
    contents = await file.read()
    print(contents)

    if filename is None:
        filename = file.filename
    if isinstance(description, str):
        file_type = filename.split('.')[-1]  # getting extension
        if file_type in ['csv', 'parquet']:
            data_type = 'pandas.DataFrame'
        elif file_type in [
                'pdb', 'sdf', 'fasta', 'fastq', 'sdf', 'txt', 'xml',
                'pdbqt', 'smi', 'smiles', 'cxsmiles', 'json'
        ]:
            data_type = 'text/plain'
        elif file_type in ['dcd', 'bz2', 'zip', 'onnx', 'hdf5']:
            data_type = 'binary'
        elif file_type in ['png']:
            data_type = 'png'
        else:
            data_type = ''
        card: Union[DataCard, dict] = DataCard(address='',
                                               file_type=file_type,
                                               data_type=data_type,
                                               description=description)
    else:
        card = description

    address = _upload_data(profile_name, project_name, filename, contents, card, backend=backend)
    return {"dataset_address": address}
