"""
V1 API Router - Aggregates all API endpoints under /v1 prefix.

This module creates a versioned API router that includes all sub-routers:
- data: Data upload/download endpoints
- primitives: ML primitive endpoints (featurize, train, evaluate, infer)
- jobs: Job status and management endpoints
"""
from fastapi import APIRouter

from deepchem_server.api.routers.data import router as data_router
from deepchem_server.api.routers.jobs import router as jobs_router
from deepchem_server.api.routers.primitives import router as primitives_router


router = APIRouter(prefix="/v1")

router.include_router(data_router)
router.include_router(primitives_router)
router.include_router(jobs_router)
