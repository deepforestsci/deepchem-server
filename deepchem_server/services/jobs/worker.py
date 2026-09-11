"""
RQ Worker function for executing jobs.

This module contains the function that RQ workers call to execute jobs.
It handles:
- Loading job from database
- Updating status to RUNNING
- Executing via ComputeWorkflow
- Updating to COMPLETED/FAILED
- Retry logic
"""
import json
import logging
import os
from datetime import datetime
from typing import Any

from deepchem_server.core.common import config
from deepchem_server.core.primitives.compute import ComputeWorkflow
from deepchem_server.services.jobs.models import Job, JobStatus, get_session
from deepchem_server.services.jobs.utils import init_datastore


logger = logging.getLogger(__name__)


def _setup_datastore(profile: str, project: str) -> None:
    """Set up the DeepchemDatastore for this worker.
    
    Requires DATASTORE_URL environment variable to be set.
    
    Parameters
    ----------
    profile : str
        Profile name
    project : str
        Project name
        
    Raises
    ------
    RuntimeError
        If DATASTORE_URL is not set
    """
    datastore = init_datastore(profile_name=profile, project_name=project)
    config.set_datastore(datastore)
    logger.info(f"Using datastore service: {os.getenv('DATASTORE_URL')}")


def execute_job(job_id: str) -> Any:
    """Execute a job by ID.
    
    This function is called by RQ workers.
    
    Parameters
    ----------
    job_id : str
        Job ID to execute
        
    Returns
    -------
    Any
        Result from ComputeWorkflow.execute()
        
    Raises
    ------
    Exception
        If execution fails
    """
    logger.info(f"Starting job execution: {job_id}")

    session = get_session()
    db_job = None
    try:
        # Load job from database
        db_job = session.query(Job).filter(Job.id == job_id).first()
        if not db_job:
            raise ValueError(f"Job not found: {job_id}")

        # Update status to RUNNING
        db_job.status = JobStatus.RUNNING
        db_job.started_at = datetime.utcnow()
        session.commit()

        # Parse program
        program = json.loads(db_job.program_json)
        profile = db_job.profile
        project = db_job.project

        logger.info(f"Executing program: {db_job.program_name}")

        # Setup datastore
        _setup_datastore(profile, project)

        # Execute workflow
        workflow = ComputeWorkflow(program)
        result = workflow.execute()

        # Update to COMPLETED
        db_job.status = JobStatus.COMPLETED
        db_job.result = json.dumps(str(result)) if result else None
        db_job.completed_at = datetime.utcnow()
        session.commit()

        logger.info(f"Job completed successfully: {job_id}")
        return result

    except Exception as e:
        logger.error(f"Job failed: {job_id} - {e}")

        # Update job with error
        if db_job:
            db_job.retry_count += 1

            if db_job.retry_count >= db_job.max_retries:
                db_job.status = JobStatus.FAILED
                db_job.error = str(e)
                db_job.completed_at = datetime.utcnow()
                logger.info(f"Job permanently failed after {db_job.retry_count} attempts: {job_id}")
            else:
                # Will be retried by RQ
                db_job.status = JobStatus.QUEUED
                db_job.error = f"Attempt {db_job.retry_count}: {str(e)}"
                logger.info(f"Job will retry (attempt {db_job.retry_count}/{db_job.max_retries}): {job_id}")

            session.commit()

        raise

    finally:
        session.close()
