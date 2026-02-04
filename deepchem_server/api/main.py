"""
DeepChem Server API Application.

Entry point for the FastAPI server. This module creates the main FastAPI app
and includes all versioned API routers.

Usage:
    uvicorn deepchem_server.api.main:app --host 0.0.0.0 --port 8000
"""
import logging

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from deepchem_server.api.routers import v1_router


logger = logging.getLogger("deepchem_server")
logger.setLevel(logging.INFO)

app = FastAPI(
    title="DeepChem Server",
    description="REST API for DeepChem computational chemistry workflows",
    version="1.0.0",
)

# Include versioned router
app.include_router(v1_router)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.on_event("startup")
async def on_startup():
    """Startup event handler."""
    pass


@app.get("/healthcheck")
async def perform_healthcheck():
    """
    HealthCheck endpoint to check server status.
    
    Returns
    -------
    JSONResponse
        Status OK response
    """
    return JSONResponse(status_code=200, content={"status": "ok"})
