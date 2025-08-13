import logging

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from deepchem_server.routers import data, primitives

logger = logging.getLogger("backend_logs")
logger.setLevel(logging.INFO)

app = FastAPI()

app.include_router(primitives.router)
app.include_router(data.router)

app.add_middleware(
    CORSMiddleware,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.on_event("startup")
async def on_startup():
    pass


@app.get("/healthcheck")
async def perform_healthcheck():
    """
    HealthCheck endpoint to check server status
    """
    return JSONResponse(status_code=200, content={"status": "ok"})
