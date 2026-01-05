"""
Jobs Router - Endpoints for job status and management.
"""
import logging
import os
from typing import Any, Dict, Optional

from fastapi import APIRouter, Depends, HTTPException, Query

from deepchem_server.services.jobs import JobQueue, JobStatus


logger = logging.getLogger(__name__)

# Singleton queue instance
_job_queue: Optional[JobQueue] = None


def get_job_queue() -> JobQueue:
    """Get or create JobQueue singleton."""
    global _job_queue
    if _job_queue is None:
        redis_url = os.getenv("REDIS_URL", "redis://localhost:6379/0")
        database_url = os.getenv("JOBS_DATABASE_URL")
        _job_queue = JobQueue(redis_url=redis_url, database_url=database_url)
    return _job_queue


# API Router
router = APIRouter(prefix="/jobs", tags=["jobs"])


@router.get("/{job_id}")
async def get_job(
        job_id: str,
        queue: JobQueue = Depends(get_job_queue),
) -> Dict[str, Any]:
    """Get job details by ID.
    
    Parameters
    ----------
    job_id : str
        Job ID
        
    Returns
    -------
    dict
        Job details including status, result, errors
    """
    job = queue.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")
    return job


@router.get("")
async def list_jobs(
        profile: Optional[str] = Query(None),
        project: Optional[str] = Query(None),
        status: Optional[str] = Query(None),
        limit: int = Query(100, le=1000),
        queue: JobQueue = Depends(get_job_queue),
) -> Dict[str, Any]:
    """List jobs with optional filters.
    
    Parameters
    ----------
    profile : str, optional
        Filter by profile
    project : str, optional
        Filter by project
    status : str, optional
        Filter by status (pending, queued, running, completed, failed)
    limit : int
        Maximum results (default 100, max 1000)
        
    Returns
    -------
    dict
        List of jobs
    """
    status_enum = None
    if status:
        try:
            status_enum = JobStatus(status)
        except ValueError:
            raise HTTPException(status_code=400,
                                detail=f"Invalid status: {status}. Valid: pending, queued, running, completed, failed")

    jobs = queue.get_jobs(
        profile=profile,
        project=project,
        status=status_enum,
        limit=limit,
    )
    return {"jobs": jobs, "count": len(jobs)}


@router.delete("/{job_id}")
async def cancel_job(
        job_id: str,
        queue: JobQueue = Depends(get_job_queue),
) -> Dict[str, Any]:
    """Cancel a pending or queued job.
    
    Parameters
    ----------
    job_id : str
        Job ID
        
    Returns
    -------
    dict
        Cancellation result
    """
    success = queue.cancel_job(job_id)
    if not success:
        raise HTTPException(status_code=400,
                            detail=f"Cannot cancel job {job_id}. May not exist or already running/completed.")
    return {"status": "cancelled", "job_id": job_id}


@router.post("/recover")
async def recover_stale_jobs(queue: JobQueue = Depends(get_job_queue),) -> Dict[str, Any]:
    """Recover jobs that were running when workers crashed.
    
    Returns
    -------
    dict
        Number of recovered jobs
    """
    recovered = queue.recover_stale_jobs()
    return {"status": "success", "recovered_count": recovered}
