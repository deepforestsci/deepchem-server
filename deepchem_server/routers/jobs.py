from fastapi import APIRouter, HTTPException

from deepchem_server.aws.batch import get_job_status


router = APIRouter(prefix="/jobs", tags=["jobs"])


@router.get("/{job_id}")
async def job_status(job_id: str) -> dict:
    """Poll the status of an AWS Batch job submitted by a primitive endpoint."""
    result = get_job_status(job_id)
    if result["status"] == "NOT_FOUND":
        raise HTTPException(status_code=404, detail=f"Job {job_id!r} not found")
    return result
