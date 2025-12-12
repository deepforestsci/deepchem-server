"""
FastAPI REST API for DeepchemDatastore service.

Provides HTTP endpoints for CRUD operations on the datastore.
"""

import logging
import os
from typing import Any, Dict, Optional

from fastapi import APIRouter, Depends, FastAPI, File, Form, HTTPException, Query, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import Response

from deepchem_server.services.datastore.auth import verify_api_key
from deepchem_server.services.datastore.service import DatastoreService


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
        _: str = Depends(verify_api_key),
        service: DatastoreService = Depends(get_datastore_service),
) -> Dict[str, Any]:
    """List all data keys in a profile/project.

    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name

    Returns
    -------
    dict
        Response with list of keys
    """
    try:
        keys = service.list_data(profile=profile, project=project)
        return {"status": "success", "keys": keys, "count": len(keys)}
    except Exception as e:
        logger.exception(f"Error listing data: {e}")
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


# Default app instance
app = create_app()
