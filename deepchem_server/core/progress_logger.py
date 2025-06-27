import sys
import logging

logging.basicConfig(stream=sys.stdout, level=logging.INFO, force=True)
logger = logging.getLogger("progress_logger")


def log_progress(job_type: str, progress: int, message: str):
    """
    Logs progress of a job

    Parameters
    ----------
    job_type: str
        Type of job
    progress: int
        Progress of the job
    message: str
        Message to be logged

    Returns
    -------
    None
    """
    assert 0 <= progress <= 100
    logger.info(f"{job_type}: {progress}% - {message}")
