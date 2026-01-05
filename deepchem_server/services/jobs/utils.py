"""
Jobs service utilities.

Contains utilities for job execution in workers.
"""
import logging
import os
from typing import Dict

from deepchem_server.core.common import config
from deepchem_server.services.datastore.client import DatastoreClient, DeepchemDatastore


logger = logging.getLogger(__name__)

DATA_DIR = os.getenv("DATADIR", "./data")


def init_datastore(profile_name: str, project_name: str) -> DeepchemDatastore:
    """
    Initialize the DeepchemDatastore for a job.
    
    Requires DATASTORE_URL environment variable to be set.

    Parameters
    ----------
    profile_name: str
        Name of the Profile where the job is run
    project_name: str
        Name of the Project where the job is run
        
    Raises
    ------
    RuntimeError
        If DATASTORE_URL is not set
        
    Returns
    -------
    DeepchemDatastore
        Initialized datastore instance
    """
    datastore_url = os.getenv("DATASTORE_URL")
    if not datastore_url:
        raise RuntimeError("DATASTORE_URL environment variable must be set. "
                           "All operations require a running datastore service.")

    client = DatastoreClient(url=datastore_url)
    return DeepchemDatastore(
        client=client,
        profile=profile_name,
        project=project_name,
        temp_dir=DATA_DIR,
    )


def run_job_sync(profile_name: str, project_name: str, program: Dict, backend: str = 'local'):
    """
    Run a job synchronously (for workers or testing).

    Parameters
    ----------
    profile_name: str
        Name of the Profile where the job is run
    project_name: str
        Name of the Project where the job is run
    program: Dict
        Program dictionary containing program name and kwargs
    backend: str
        Backend to be used to run the job (Default: local)
        
    Returns
    -------
    Any
        Result from ComputeWorkflow.execute()
    """
    # Lazy import to avoid loading deepchem in gateway
    from deepchem_server.core.primitives.compute import ComputeWorkflow

    logger.info("beginning")
    datastore = init_datastore(profile_name=profile_name, project_name=project_name)
    config.set_datastore(datastore)
    workflow = ComputeWorkflow(program)
    try:
        output = workflow.execute()
    except Exception as e:
        logger.error(f"Error executing workflow: {e}")
        raise e
    return output
