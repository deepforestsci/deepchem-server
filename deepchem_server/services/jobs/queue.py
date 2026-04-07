"""
Redis Queue wrapper for job management.

Handles job submission, status tracking, and database synchronization.
"""
import json
import logging
import os
import uuid
from datetime import datetime
from typing import Any, Dict, List, Optional

from redis import Redis
from rq import Queue, Retry
from rq.job import Job as RQJob

from deepchem_server.services.jobs.models import Job, JobStatus, get_session, init_db


logger = logging.getLogger(__name__)

# Redis configuration
REDIS_URL = os.getenv("REDIS_URL", "redis://localhost:6379/0")
DEFAULT_QUEUE_NAME = "deepchem_jobs"
DEFAULT_TIMEOUT = 3600  # 1 hour


class JobQueue:
    """Redis Queue manager with persistent job state tracking.
    
    Provides job submission, status tracking, and recovery mechanisms.
    
    Parameters
    ----------
    redis_url : str
        Redis connection URL
    queue_name : str
        Name of the RQ queue
    database_url : str, optional
        MySQL database URL for job persistence
    """

    _job_queue_name: str = "deepchem_server.services.jobs.worker.execute_job"
    _result_ttl: int = 86400 * 7  # 7 days
    _failure_ttl: int = 86400 * 7

    def __init__(
        self,
        redis_url: str = REDIS_URL,
        queue_name: str = DEFAULT_QUEUE_NAME,
        database_url: Optional[str] = None,
    ) -> None:
        self.redis_conn = Redis.from_url(redis_url)
        self.queue = Queue(queue_name, connection=self.redis_conn)
        init_db(database_url)
        logger.info(f"JobQueue initialized: queue={queue_name}")

    def enqueue(
        self,
        profile: str,
        project: str,
        program: Dict[str, Any],
        max_retries: int = 3,
        timeout: int = DEFAULT_TIMEOUT,
    ) -> str:
        """Submit a job to the queue.
        
        Parameters
        ----------
        profile : str
            Profile name for datastore
        project : str
            Project name for datastore
        program : dict
            Program configuration (must contain 'program_name')
        max_retries : int
            Maximum retry attempts
        timeout : int
            Job timeout in seconds
            
        Returns
        -------
        str
            Job ID
        """
        job_id = str(uuid.uuid4())
        program_name = program.get("program_name", "unknown")

        # Create database record
        session = get_session()
        try:
            db_job = Job(
                id=job_id,
                program_name=program_name,
                program_json=json.dumps(program),
                profile=profile,
                project=project,
                max_retries=max_retries,
            )
            session.add(db_job)
            session.commit()
            logger.info(f"Created job record: {job_id}")
        except Exception as e:
            session.rollback()
            logger.error(f"Failed to create job record: {e}")
            raise
        finally:
            session.close()

        # Enqueue to Redis
        try:
            _ = self.queue.enqueue(self._job_queue_name,
                                   job_id,
                                   job_timeout=timeout,
                                   job_id=job_id,
                                   result_ttl=self._result_ttl,
                                   failure_ttl=self._failure_ttl,
                                   retry=Retry(max=max_retries))

            # Update status to QUEUED
            self._update_job_status(job_id, JobStatus.QUEUED)
            logger.info(f"Enqueued job: {job_id}")

        except Exception as e:
            self._update_job_status(job_id, JobStatus.FAILED, error=str(e))
            logger.error(f"Failed to enqueue job: {e}")
            raise

        return job_id

    def get_job(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Get job details by ID.
        
        Parameters
        ----------
        job_id : str
            Job ID
            
        Returns
        -------
        dict or None
            Job details or None if not found
        """
        session = get_session()
        try:
            db_job = session.query(Job).filter(Job.id == job_id).first()
            if db_job:
                return db_job.to_dict()
            return None
        finally:
            session.close()

    def get_jobs(
        self,
        profile: Optional[str] = None,
        project: Optional[str] = None,
        status: Optional[JobStatus] = None,
        limit: int = 100,
    ) -> List[Dict[str, Any]]:
        """List jobs with optional filters.
        
        Parameters
        ----------
        profile : str, optional
            Filter by profile
        project : str, optional
            Filter by project
        status : JobStatus, optional
            Filter by status
        limit : int
            Maximum number of results
            
        Returns
        -------
        list of dict
            List of job details
        """
        session = get_session()
        try:
            query = session.query(Job)

            if profile:
                query = query.filter(Job.profile == profile)
            if project:
                query = query.filter(Job.project == project)
            if status:
                query = query.filter(Job.status == status)

            jobs = query.order_by(Job.created_at.desc()).limit(limit).all()
            return [job.to_dict() for job in jobs]
        finally:
            session.close()

    def cancel_job(self, job_id: str) -> bool:
        """Cancel a pending or queued job.
        
        Parameters
        ----------
        job_id : str
            Job ID
            
        Returns
        -------
        bool
            True if cancelled successfully
        """
        session = get_session()
        try:
            db_job = session.query(Job).filter(Job.id == job_id).first()
            if not db_job:
                return False

            if db_job.status not in (JobStatus.PENDING, JobStatus.QUEUED):
                return False

            # Remove from Redis queue
            try:
                rq_job = RQJob.fetch(job_id, connection=self.redis_conn)
                rq_job.cancel()
            except Exception:
                pass  # Job may not be in queue

            db_job.status = JobStatus.FAILED
            db_job.error = "Cancelled by user"
            db_job.completed_at = datetime.utcnow()
            session.commit()
            return True
        finally:
            session.close()

    def recover_stale_jobs(self) -> int:
        """Recover jobs that were running when workers crashed.
        
        Finds jobs marked as RUNNING and re-queues them.
        
        Returns
        -------
        int
            Number of jobs recovered
        """
        session = get_session()
        try:
            stale_jobs = session.query(Job).filter((Job.status == JobStatus.RUNNING) or
                                                   (Job.status == JobStatus.QUEUED)).all()
            recovered = 0

            for job in stale_jobs:
                if job.retry_count < job.max_retries:
                    job.retry_count += 1
                    job.status = JobStatus.QUEUED

                    # Re-enqueue to Redis
                    self.queue.enqueue(self._job_queue_name,
                                       job.id,
                                       job_id=job.id,
                                       job_timeout=DEFAULT_TIMEOUT,
                                       result_ttl=self._result_ttl,
                                       failure_ttl=self._failure_ttl,
                                       retry=Retry(max=job.max_retries))
                    recovered += 1
                    logger.info(f"Recovered stale job: {job.id} (attempt {job.retry_count})")
                else:
                    job.status = JobStatus.FAILED
                    job.error = "Max retries exceeded after worker crash"
                    job.completed_at = datetime.utcnow()

            session.commit()
            logger.info(f"Recovered {recovered} stale jobs")
            return recovered
        finally:
            session.close()

    def _update_job_status(
        self,
        job_id: str,
        status: JobStatus,
        result: Optional[str] = None,
        error: Optional[str] = None,
    ) -> None:
        """Update job status in database."""
        session = get_session()
        try:
            job = session.query(Job).filter(Job.id == job_id).first()
            if job:
                job.status = status
                if status == JobStatus.RUNNING:
                    job.started_at = datetime.utcnow()
                if status in (JobStatus.COMPLETED, JobStatus.FAILED):
                    job.completed_at = datetime.utcnow()
                if result:
                    job.result = result
                if error:
                    job.error = error
                session.commit()
        finally:
            session.close()
