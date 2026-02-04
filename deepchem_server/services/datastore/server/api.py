"""
FastAPI REST API for DeepchemDatastore service.

Provides HTTP endpoints for CRUD operations on the datastore.
"""
import json
import logging
import os
from typing import Annotated, Any, Dict, List, Optional

from fastapi import (
    APIRouter,
    Depends,
    FastAPI,
    File,
    Form,
    HTTPException,
    Query,
    Body,
    UploadFile,
)
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import Response

from deepchem_server.services.datastore.server.auth import verify_api_key
from deepchem_server.services.datastore.server.service import DatastoreService


logger = logging.getLogger(__name__)

DEFAULT_BASE_DIR = os.getenv("DATASTORE_BASE_DIR", "./datastore_data")
DEFAULT_PORT = int(os.getenv("DATASTORE_PORT", "8081"))

_datastore_service: Optional[DatastoreService] = None


def get_datastore_service() -> DatastoreService:
    """Get or create the datastore service singleton."""
    global _datastore_service
    if _datastore_service is None:
        base_dir = os.getenv("DATASTORE_BASE_DIR", DEFAULT_BASE_DIR)
        _datastore_service = DatastoreService(base_dir)
    return _datastore_service


def create_app(base_dir: Optional[str] = None) -> FastAPI:
    """Create the FastAPI application.

    Parameters
    ----------
    base_dir : str, optional
        Base directory for storage. Defaults to DATASTORE_BASE_DIR env var.

    Returns
    -------
    FastAPI
        The configured FastAPI application
    """
    global _datastore_service

    if base_dir:
        _datastore_service = DatastoreService(base_dir)

    app = FastAPI(
        title="DeepchemDatastore Service",
        description="Standalone file storage service following DeepchemAddress pattern",
        version="1.0.0",
    )

    # CORS middleware for development
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # Include routes
    app.include_router(router)

    return app


router = APIRouter(prefix="/api/v1", tags=["datastore"])


@router.get("/healthcheck")
async def healthcheck() -> Dict[str, str]:
    """Health check endpoint (no auth required)."""
    return {"status": "ok", "service": "datastore"}


@router.post("/data/{profile}/{project}/{key:path}")
async def upload_data(
        profile: str,
        project: str,
        key: str,
        file: UploadFile = File(...),
        card: Optional[str] = Form(None),
        kind: str = Form("data"),
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Dict[str, Any]:
    """Upload data to the datastore.

    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name
    key : str
        File key/path
    file : UploadFile
        File to upload
    card : str, optional
        JSON string of metadata card
    kind : str
        Type of object ('data' or 'model')

    Returns
    -------
    dict
        Response with address and status
    """
    try:
        data = await file.read()
        card_dict = None
        if card:
            import json

            card_dict = json.loads(card)

        address = service.upload_data(
            profile=profile,
            project=project,
            key=key,
            data=data,
            card=card_dict,
            kind=kind,
        )

        return {
            "status": "success",
            "address": address,
            "size": len(data),
        }
    except Exception as e:
        logger.exception(f"Error uploading data: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/dir/{profile}/{project}/{dirname}")
async def create_directory(
        profile: str,
        project: str,
        dirname: str,
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Dict[str, Any]:
    """Create a directory in the datastore."""
    try:
        address = service.create_directory(
            profile=profile,
            project=project,
            dirname=dirname,
        )
        return {"status": "success", "address": address}
    except Exception as e:
        logger.exception(f"Error creating directory: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/data/{profile}/{project}/{key:path}")
async def get_data(
        profile: str,
        project: str,
        key: str,
        fetch_sample: bool = Query(False),
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Response:
    """Get data from the datastore.

    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name
    key : str
        File key/path
    fetch_sample : bool
        Whether to fetch only a sample

    Returns
    -------
    Response
        File content as bytes
    """
    try:
        data = service.get_data(
            profile=profile,
            project=project,
            key=key,
            fetch_sample=fetch_sample,
        )

        # Determine content type based on extension
        content_type = "application/octet-stream"
        if key.endswith(".csv"):
            content_type = "text/csv"
        elif key.endswith(".json"):
            content_type = "application/json"
        elif key.endswith(".txt"):
            content_type = "text/plain"

        return Response(content=data, media_type=content_type)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Object not found: {profile}/{project}/{key}")
    except Exception as e:
        logger.exception(f"Error getting data: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/move/{profile}/{project}/{key:path}")
async def move_object(
        profile: str,
        project: str,
        key: str,
        destination: Annotated[str, Body(embed=True)],
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Dict[str, Any]:
    """Move an object to a new location.

    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name
    key : str
        File key/path
    destination : str
        New location address

    Returns
    -------
    dict
        Response with status and new address
    """
    try:
        address = service.move_object(
            profile=profile,
            project=project,
            key=key,
            destination=destination,
        )
        return {"status": "success", "address": address}
    except Exception as e:
        logger.exception(f"Error moving object: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.delete("/data/{profile}/{project}/{key:path}")
async def delete_data(
        profile: str,
        project: str,
        key: str,
        kind: str = Query("data"),
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Dict[str, Any]:
    """Delete data from the datastore.

    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name
    key : str
        File key/path
    kind : str
        Type of object ('data', 'model', 'dir')

    Returns
    -------
    dict
        Response with status
    """
    try:
        service.delete_object(
            profile=profile,
            project=project,
            key=key,
            kind=kind,
        )
        return {"status": "success", "deleted": f"{profile}/{project}/{key}"}
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Object not found: {profile}/{project}/{key}")
    except Exception as e:
        logger.exception(f"Error deleting data: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/list/{profile}/{project}")
async def list_data(
        profile: str,
        project: str,
        include_card_files: bool = Query(False),
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Dict[str, Any]:
    """List all data keys in a profile/project (excludes card files).

    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name
    include_card_files : bool
        Whether to include card files (.cdc, .cmc)
    Returns
    -------
    dict
        Response with list of keys
    """
    try:
        keys = service.list_data(profile=profile, project=project, include_card_files=include_card_files)
        return {"status": "success", "keys": keys, "count": len(keys)}
    except Exception as e:
        logger.exception(f"Error listing data: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/list-all/{profile}/{project}")
async def list_all_objects(
        profile: str,
        project: str,
        prefix: str = Query(""),
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Dict[str, Any]:
    """List all objects including card files (.cdc, .cmc).

    This endpoint is used by DeepchemDatastore._get_datastore_objects() for
    checkpoint scanning during multicore featurization restart.

    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name

    Returns
    -------
    dict
        Response with list of all objects including card files
    """
    try:
        objects = service.list_all_objects(profile=profile, project=project, prefix=prefix)
        return {"status": "success", "objects": objects, "count": len(objects)}
    except Exception as e:
        logger.exception(f"Error listing all objects: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/card/{profile}/{project}/{key:path}")
async def get_card(
        profile: str,
        project: str,
        key: str,
        kind: str = Query("data"),
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Dict[str, Any]:
    """Get the metadata card for an object.

    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name
    key : str
        File key/path
    kind : str
        Type of object ('data' or 'model')

    Returns
    -------
    dict
        Card metadata or empty dict
    """
    try:
        card = service.get_card(
            profile=profile,
            project=project,
            key=key,
            kind=kind,
        )
        return {"status": "success", "card": card}
    except Exception as e:
        logger.exception(f"Error getting card: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/exists/{profile}/{project}/{key:path}")
async def check_exists(
        profile: str,
        project: str,
        key: str,
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Dict[str, Any]:
    """Check if an object exists.

    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name
    key : str
        File key/path

    Returns
    -------
    dict
        Response with exists boolean
    """
    exists = service.exists(profile=profile, project=project, key=key)
    return {"status": "success", "exists": exists}


@router.get("/size/{profile}/{project}/{key:path}")
async def get_size(
        profile: str,
        project: str,
        key: str,
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Dict[str, Any]:
    """Get the size of an object.

    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name
    key : str
        File key/path

    Returns
    -------
    dict
        Response with size in bytes
    """
    try:
        size = service.get_file_size(profile=profile, project=project, key=key)
        return {"status": "success", "size": size}
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Object not found: {profile}/{project}/{key}")
    except Exception as e:
        logger.exception(f"Error getting size: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/directory/{profile}/{project}/{key:path}")
async def upload_directory(
        profile: str,
        project: str,
        key: str,
        files: List[UploadFile] = File(...),
        card: Optional[str] = Form(None),
        kind: str = Form("data"),
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Dict[str, Any]:
    """Upload multiple files as a directory structure.

    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name
    key : str
        Directory key/path
    files : List[UploadFile]
        Files to upload. Filenames should include relative paths.
    card : str, optional
        JSON string of metadata card
    kind : str
        Type of object ('data' or 'model')

    Returns
    -------
    dict
        Response with address and file count
    """
    try:
        file_tuples = []
        for idx, f in enumerate(files):
            data = await f.read()
            # Use filename as relative path (client sends relative paths as filenames)
            relative_path = f.filename or f"file_{idx}"
            file_tuples.append((relative_path, data))

        card_dict = None
        if card:
            card_dict = json.loads(card)

        address = service.upload_directory_files(
            profile=profile,
            project=project,
            key=key,
            files=file_tuples,
            card=card_dict,
            kind=kind,
        )

        return {
            "status": "success",
            "address": address,
            "file_count": len(file_tuples),
        }
    except Exception as e:
        logger.exception(f"Error uploading directory: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/directory-contents/{profile}/{project}/{key:path}")
async def get_directory_contents(
        profile: str,
        project: str,
        key: str,
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Dict[str, Any]:
    """List all files in a directory.

    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name
    key : str
        Directory key/path

    Returns
    -------
    dict
        Response with list of file paths
    """
    try:
        files = service.list_directory_contents(
            profile=profile,
            project=project,
            key=key,
        )
        return {"status": "success", "files": files, "count": len(files)}
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Directory not found: {profile}/{project}/{key}")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception(f"Error listing directory contents: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/download/{profile}/{project}/{key:path}")
async def download_object(
        profile: str,
        project: str,
        key: str,
        format: str = Query("auto", description="Download format: 'zip', 'direct', or 'auto'"),
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Response:
    """Download a file or directory.

    For directories:
    - format=zip: Returns a ZIP archive of the directory
    - format=direct: Returns JSON list of files (client downloads individually)
    - format=auto: ZIP for directories, direct for files

    For files:
    - Always returns the file content directly

    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name
    key : str
        File or directory key/path
    format : str
        Download format ('zip', 'direct', or 'auto')

    Returns
    -------
    Response
        File content or ZIP archive
    """
    from pathlib import Path

    try:
        storage_path = service._get_storage_path(profile, project, key)

        if not storage_path.exists():
            raise HTTPException(status_code=404, detail=f"Object not found: {profile}/{project}/{key}")

        # Handle single file
        if storage_path.is_file():
            data = storage_path.read_bytes()
            content_type = "application/octet-stream"
            if key.endswith(".csv"):
                content_type = "text/csv"
            elif key.endswith(".json"):
                content_type = "application/json"
            elif key.endswith(".txt"):
                content_type = "text/plain"

            filename = Path(key).name
            return Response(
                content=data,
                media_type=content_type,
                headers={"Content-Disposition": f'attachment; filename="{filename}"'},
            )

        # Handle directory
        if format == "direct" or (format == "auto" and not storage_path.is_dir()):
            # Return list of files for client to download individually
            import json
            files = service.list_directory_contents(profile, project, key)
            return Response(
                content=json.dumps({"files": files}),
                media_type="application/json",
            )

        # Default: return as ZIP
        zip_data = service.download_directory_as_zip(profile, project, key)
        filename = f"{Path(key).name}.zip"
        return Response(
            content=zip_data,
            media_type="application/zip",
            headers={"Content-Disposition": f'attachment; filename="{filename}"'},
        )

    except HTTPException:
        raise
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Object not found: {profile}/{project}/{key}")
    except Exception as e:
        logger.exception(f"Error downloading: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# Default app instance
app = create_app()
