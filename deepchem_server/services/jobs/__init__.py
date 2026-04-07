# Jobs Service Package
"""
Redis Job Queue with MySQL/SQLite persistence for job state tracking.

Components:
- models: SQLAlchemy models for job persistence
- queue: Job queue manager (requires redis)
- worker: RQ worker function (requires redis)
- utils: Utility functions for job execution

Usage:
    from deepchem_server.services.jobs import run_job, JobQueue, JobStatus
"""
import logging

from deepchem_server.services.jobs.models import Job, JobStatus
from deepchem_server.services.jobs.queue import JobQueue


logger = logging.getLogger(__name__)


def run_job(profile_name: str, project_name: str, program: dict, backend: str = 'local') -> str:
    """
    Submit a job to the job queue for async execution.

    Parameters
    ----------
    profile_name: str
        Name of the Profile where the job is run
    project_name: str
        Name of the Project where the job is run
    program: dict
        Program dictionary containing program name and kwargs
    backend: str
        Backend to be used to run the job (Default: local, ignored)
        
    Returns
    -------
    str
        Job ID for tracking the job status
    """
    queue = JobQueue()
    job_id = queue.enqueue(
        profile=profile_name,
        project=project_name,
        program=program,
    )
    logger.info(f"Job submitted: {job_id}")
    return job_id


__all__ = ["Job", "JobStatus", "JobQueue", "run_job"]
