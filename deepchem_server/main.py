from fastapi import FastAPI
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from deepchem_server.routers import primitives, data

import logging

logger = logging.getLogger("backend_logs")
logger.setLevel(logging.INFO)

app = FastAPI()

app.include_router(primitives.router)
app.include_router(data.router)

origins = [
    "http://localhost:4200",
    "https://localhost:4200",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.on_event("startup")
async def on_startup():
    pass


@app.get("/healthcheck")
async def perform_healthcheck():
    return JSONResponse(status_code=200, content={"status": "ok"})
