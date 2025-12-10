import json
import os

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse

from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.common.datastore import DiskDataStore


app = FastAPI(title="DeepChem Datastore API", version="1.0.0")

# Config - this could be moved to environment variables or a config file
DATA_DIR = os.getenv("DATADIR", "./data")


def get_datastore(profile_name: str, project_name: str) -> DiskDataStore:
    """Get or create a DiskDataStore instance for the given profile/project."""
    return DiskDataStore(profile_name=profile_name, project_name=project_name, basedir=DATA_DIR)


def _get_data_type_from_file_type(file_type: str) -> str:
    """Map file type to appropriate data type."""
    type_mapping = {
        "csv": "pandas.DataFrame",
        "json": "json",
        "pdb": "text/plain",
        "sdf": "text/plain",
        "fasta": "text/plain",
        "fastq": "text/plain",
        "txt": "text/plain",
        "xml": "text/plain",
        "pdbqt": "text/plain",
        "smi": "text/plain",
        "smiles": "text/plain",
        "cxsmiles": "text/plain",
        "png": "png",
        "dcd": "binary",
        "bz2": "binary",
        "zip": "binary",
        "onnx": "binary",
        "hdf5": "binary",
        "log": "text/plain",
    }
    return type_mapping.get(file_type, "binary")


@app.post("/v1/files/{profile_name}/{project_name}")
async def upload_file(
    profile_name: str,
    project_name: str,
    file: UploadFile = File(...),
):
    """Upload a file to a specific project."""
    datastore = get_datastore(profile_name, project_name)

    file_type = file.filename.split(".")[-1] if "." in file.filename else "txt"
    data_type = _get_data_type_from_file_type(file_type)

    card = DataCard(
        address="",
        file_type=file_type,
        data_type=data_type,
        description=f"Uploaded file {file.filename}",
    )

    content = await file.read()

    try:
        dataset_address = datastore.upload_data(
            datastore_filename=file.filename, filename=content, card=card
        )
        return JSONResponse(status_code=201, content={"address": dataset_address})
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/v1/files/{profile_name}/{project_name}/{filename:path}")
async def download_file(profile_name: str, project_name: str, filename: str):
    """Download a file."""
    datastore = get_datastore(profile_name, project_name)
    address = f"deepchem://{profile_name}/{project_name}/{filename}"

    if not datastore.exists(address):
        raise HTTPException(status_code=404, detail="File not found")

    file_path = os.path.join(datastore.storage_loc, filename)
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="File not found on disk")

    return FileResponse(file_path)


@app.head("/v1/files/{profile_name}/{project_name}/{filename:path}")
async def check_file_exists(profile_name: str, project_name: str, filename: str):
    """Check if a file exists."""
    datastore = get_datastore(profile_name, project_name)
    address = f"deepchem://{profile_name}/{project_name}/{filename}"

    if datastore.exists(address):
        return JSONResponse(status_code=200, content={"exists": True})
    else:
        raise HTTPException(status_code=404, detail="File not found")


@app.delete("/v1/files/{profile_name}/{project_name}/{filename:path}")
async def delete_file(profile_name: str, project_name: str, filename: str):
    """Delete a file."""
    datastore = get_datastore(profile_name, project_name)
    address = f"deepchem://{profile_name}/{project_name}/{filename}"

    if not datastore.exists(address):
        raise HTTPException(status_code=404, detail="File not found")

    try:
        datastore.delete_object(address)
        return JSONResponse(status_code=200, content={"deleted": True})
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/v1/cards/{profile_name}/{project_name}/{filename:path}")
async def get_card(profile_name: str, project_name: str, filename: str):
    """Get card metadata for a file."""
    datastore = get_datastore(profile_name, project_name)
    address = f"deepchem://{profile_name}/{project_name}/{filename}"

    card = datastore.get_card(address)
    if card:
        return JSONResponse(content=json.loads(card.to_json()))
    else:
        raise HTTPException(status_code=404, detail="Card not found")
