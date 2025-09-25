import logging
import os
from typing import Dict

from deepchem_server.core import config
from deepchem_server.core.compute import ComputeWorkflow
from deepchem_server.core.datastore import DataStore, DiskDataStore


logger = logging.getLogger(__name__)

DATA_DIR = os.getenv("DATADIR", "./data")


def _init_datastore(profile_name: str, project_name: str, backend='local') -> DataStore:
    """
    Function to initialise the datastore in DATA_DIR

    Parameters
    ----------
    profile_name: str
        Name of the Profile where the job is run
    project_name: str
        Name of the Project where the job is run
    backend: str
        Backend to be used to run the job (Default: local)
    """
    if backend == 'local':
        datastore: DataStore = DiskDataStore(profile_name=profile_name, project_name=project_name, basedir=DATA_DIR)
    else:
        raise NotImplementedError(f"{backend} backend not implemented")
    return datastore


def run_job(profile_name: str, project_name: str, program: Dict, backend: str = 'local'):
    """
    Function to run jobs based on the submitted program

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
    """
    if backend == 'local':
        logger.info("beginning")
        datastore: DataStore = _init_datastore(profile_name=profile_name, project_name=project_name, backend=backend)
        config.set_datastore(datastore)  # type: ignore
        workflow = ComputeWorkflow(program)
        try:
            output = workflow.execute()
        except Exception as e:
            logger.error(f"Error executing workflow: {e}")
            raise e
        return output
    else:
        raise NotImplementedError(f"{backend} backend not implemented")


def _upload_data(profile_name, project_name, datastore_filename, contents, data_card, backend='local'):
    """
    A wrapper method to the server for creating DataStore object and using
    it to upload data files

    Parameters
    ----------
    profile_name: str
        Name of the Profile where the job is run
    project_name: str
        Name of the Project where the job is run
    datastore_filename: str
        The file name in which data is to be written in DataStore
    contents: object
        The filepath in disk or file object in memory from which
        data will be read for writing to datastore
    data_card: DataCard
        data card for the file
    """
    datastore = _init_datastore(profile_name=profile_name, project_name=project_name, backend=backend)
    import tempfile
    tempdir = tempfile.TemporaryDirectory()
    temppath = os.path.join(tempdir.name, datastore_filename.replace('/', '_'))
    with open(temppath, 'wb') as f:
        f.write(contents)
    dataset_address = datastore.upload_data(datastore_filename=datastore_filename, filename=temppath, card=data_card)
    return dataset_address


def parse_boolean_none_values_from_kwargs(kwargs: Dict) -> Dict:
    """
    Parse boolean values from kwargs and convert 'None' to None.

    Parameters
    ----------
    kwargs : Dict
        Dictionary of values to be parsed.

    Returns
    -------
    Dict
        Dictionary with boolean values and None where applicable.
    """
    parsed_kwargs: Dict = {}
    for key, value in kwargs.items():
        if isinstance(value, str):
            if value.lower() == "true":
                parsed_kwargs[key] = True
            elif value.lower() == "false":
                parsed_kwargs[key] = False
            elif value.lower() == "none":
                parsed_kwargs[key] = None
            else:
                parsed_kwargs[key] = value
        else:
            parsed_kwargs[key] = value
    return parsed_kwargs
