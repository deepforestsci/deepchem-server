import os
from deepchem_server.core import config
from deepchem_server.core.datastore import DiskDataStore
from deepchem_server.core.compute import ComputeWorkflow


def _init_datastore(profile_name: str,
                    project_name: str,
                    backend='local'):
    if backend == 'local':
        datastore = DiskDataStore(profile_name=profile_name,
                                  project_name=project_name)
    else:
        raise NotImplementedError(f"{backend} backend not implemented")
    return datastore


def run_job(profile_name: str,
            project_name: str,
            program: str,
            backend: str = 'local'):
    """
    """
    if backend == 'local':
        print("beginning")
        datastore = _init_datastore(profile_name=profile_name,
                                    project_name=project_name,
                                    backend=backend)
        config.set_datastore(datastore)
        workflow = ComputeWorkflow(program)
        try:
            output = workflow.execute()
        except Exception as e:
            print(e)
            output = e
        return output
    else:
        raise NotImplementedError(f"{backend} backend not implemented")


def _upload_data(profile_name, project_name, datastore_filename, contents, data_card, backend='local'):
    """
  A wrapper method to the server for creating S3DataStore object and using
  it to upload data files

  Parameters
  ----------
  project: Project
  backend: str
    user backend ('aws', 'local')
  working_dir: str
    preferred working_dir of the user
  datastore_filename: str
    The file name in which data is to be written in S3DataStore
  contents: object
    The filepath in disk or file object in memory from which
  data will be read for writing to datastore
  data_card: dict
    data card for the file
  """
    datastore = _init_datastore(profile_name=profile_name,
                                    project_name=project_name,
                                    backend=backend)
    import tempfile
    tempdir = tempfile.TemporaryDirectory()
    temppath = os.path.join(tempdir.name, datastore_filename.replace('/', '_'))
    with open(temppath, 'wb') as f:
        f.write(contents)
    dataset_address = datastore.upload_data(
        datastore_filename=datastore_filename,
        filename=temppath,
        card=data_card)
    return dataset_address
