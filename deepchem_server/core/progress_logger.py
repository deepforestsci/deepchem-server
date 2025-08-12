import logging
import sys

logging.basicConfig(stream=sys.stdout, level=logging.INFO, force=True)
logger = logging.getLogger("progress_logger")

def log_progress(job_type: str, progress: int, message: str) -> None:
    """Log progress of a job.

    Parameters
    ----------
    job_type : str
        Type of job being tracked.
    progress : int
        Progress of the job as a percentage (0-100).
    message : str
        Message to be logged along with the progress.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If progress is not between 0 and 100 inclusive.
    """
    assert 0 <= progress <= 100
    logger.info(f"{job_type}: {progress}% - {message}")
