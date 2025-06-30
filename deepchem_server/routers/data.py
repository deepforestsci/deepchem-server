import logging
from typing import Union, Dict, List, Optional, Any
from fastapi import APIRouter, UploadFile, File, Form, Query, HTTPException
from deepchem_server.core.cards import DataCard
from deepchem_server.core.datastore import DiskDataStore
from deepchem_server.utils import _upload_data, _init_datastore
from thefuzz import process

logger = logging.getLogger("backend_logs")
logger.setLevel(logging.INFO)

router = APIRouter(prefix="/data", tags=["data"])


@router.post("/uploaddata")
async def upload_data(file: UploadFile = File(...),
                      profile_name: str = Form(None),
                      project_name: str = Form(None),
                      filename: str = Form(None),
                      description: Union[Dict, str] = Form(None),
                      backend='local') -> Dict:
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
    print(contents)

    if filename is None:
        filename = file.filename
    if isinstance(description, str):
        file_type = filename.split('.')[-1]  # getting extension
        if file_type in ['csv', 'parquet']:
            data_type = 'pandas.DataFrame'
        elif file_type in [
                'pdb', 'sdf', 'fasta', 'fastq', 'sdf', 'txt', 'xml', 'pdbqt',
                'smi', 'smiles', 'cxsmiles', 'json'
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
    if address is None:
        raise HTTPException(status_code=500, detail="Failed to upload data")

    return {"dataset_address": address}


@router.get("/list")
async def list_datastore(
    profile_name: str = Query(..., description="Name of the Profile"),
    project_name: str = Query(..., description="Name of the Project"),
    backend: str = Query("local", description="Backend to be used"),
) -> Dict[str, Any]:
    """
    List all data objects in the datastore

    Parameters
    ----------
    profile_name: str
        Name of the Profile where the data is stored
    project_name: str
        Name of the Project where the data is stored
    backend: str
        Backend to be used (Default: local)

    Returns
    -------
    Dict with list of datastore objects
    """
    try:
        datastore = _init_datastore(
            profile_name=profile_name, project_name=project_name, backend=backend
        )

        if not isinstance(datastore, DiskDataStore):
            raise HTTPException(status_code=500, detail="Unsupported datastore type")

        all_objects = datastore._get_datastore_objects(datastore.storage_loc)

        datastore_objects = []
        for obj in all_objects:
            full_address = f"deepchem://{profile_name}/{project_name}/{obj}"
            datastore_objects.append(full_address)

        return {
            "objects": datastore_objects,
            "count": len(datastore_objects),
            "profile": profile_name,
            "project": project_name,
        }

    except Exception as e:
        logger.error(f"Error listing datastore: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Failed to list datastore: {str(e)}")


@router.get("/search")
async def search_datastore(
    query: str = Query(..., description="Search query"),
    profile_name: str = Query(..., description="Name of the Profile"),
    project_name: str = Query(..., description="Name of the Project"),
    limit: int = Query(10, description="Maximum number of results to return"),
    score_cutoff: int = Query(60, description="Minimum similarity score (0-100)"),
    backend: str = Query("local", description="Backend to be used"),
) -> Dict[str, Any]:
    """
    Search data objects in the datastore using fuzzy string matching

    Parameters
    ----------
    query: str
        Search query string
    profile_name: str
        Name of the Profile where the data is stored
    project_name: str
        Name of the Project where the data is stored
    limit: int
        Maximum number of results to return (Default: 10)
    score_cutoff: int
        Minimum similarity score from 0-100 (Default: 60)
    backend: str
        Backend to be used (Default: local)

    Returns
    -------
    Dict with search results including matched objects and their scores
    """
    try:
        datastore = _init_datastore(
            profile_name=profile_name, project_name=project_name, backend=backend
        )

        if not isinstance(datastore, DiskDataStore):
            raise HTTPException(status_code=500, detail="Unsupported datastore type")

        all_objects = datastore._get_datastore_objects(datastore.storage_loc)

        if not all_objects:
            return {
                "results": [],
                "query": query,
                "count": 0,
                "profile": profile_name,
                "project": project_name,
            }

        matches = process.extract(query, all_objects, limit=limit)

        search_results = []
        for match, score in matches:
            if score >= score_cutoff:
                full_address = f"deepchem://{profile_name}/{project_name}/{match}"
                search_results.append(
                    {"object_name": match, "full_address": full_address, "similarity_score": score}
                )

        return {
            "results": search_results,
            "query": query,
            "count": len(search_results),
            "profile": profile_name,
            "project": project_name,
            "score_cutoff": score_cutoff,
        }

    except Exception as e:
        logger.error(f"Error searching datastore: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Failed to search datastore: {str(e)}")
