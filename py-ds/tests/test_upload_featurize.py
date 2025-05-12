import requests
from http import HTTPStatus
from requests_toolbelt.multipart.encoder import MultipartEncoder

BACKEND_URL = "http://0.0.0.0:8000/"
_SESSION = None

def uploadfile(profile_name, project_name, datastore_filename, filename, description=''):
    global _SESSION
    m = MultipartEncoder(
        fields={
            'filename': datastore_filename,
            'profile_name': profile_name,
            'project_name': project_name,
            'file': (filename, open(filename, 'rb'), 'text/plain'),
            'description': description
        })
    response = _SESSION.post(BACKEND_URL + "data/uploaddata",
                         data=m,
                         headers={'Content-Type': m.content_type})
    if response.status_code == HTTPStatus.OK.value:
        response = response.json()
        return response['dataset_address']
    elif response.status_code == HTTPStatus.INTERNAL_SERVER_ERROR.value:
        print("Internal server error")
        return None
    elif response.status_code >= 400:
        detail = response.json()
        print(detail['detail'])
        return None
    return None


def featurize(profile_name,
              project_name,
              dataset_address,
              featurizer,
              output,
              dataset_column,
              label_column,
              feat_kwargs = {},
              ):
    global _SESSION
    params = {
        'profile_name': profile_name, 
        'project_name': project_name, 
        'dataset_address': dataset_address,
        'featurizer': featurizer,
        'output': output,
        'dataset_column': dataset_column,
        'label_column': label_column,
    }

    json_params = {'feat_kwargs': feat_kwargs}
    api_path = "primitive/featurize"
    api_endpoint = BACKEND_URL+api_path
    response = _SESSION.post(api_endpoint,
                                    params=params, json=json_params)
    if response.status_code == HTTPStatus.OK.value:
        response = response.json()
        return response['featurized_file_address']
    return None


def test_upload_csv():
    global _SESSION
    _SESSION = requests.Session()

    profile_name = 'test-profile'
    project_name = 'test-project'

    dataset_address = uploadfile(
                profile_name = profile_name,
                project_name = project_name,
                datastore_filename='zinc10_sample.csv',
                filename='./assets/zinc10.csv',
                description='Sample test csv file')

    if dataset_address is None:
        raise Exception("CSV file upload failed")

    assert dataset_address == f"deepchem://{profile_name}/{project_name}/zinc10_sample.csv"


def test_featurize():
    global _SESSION
    _SESSION = requests.Session()

    profile_name = 'test-profile'
    project_name = 'test-project'

    dataset_address = uploadfile(
                profile_name = profile_name,
                project_name = project_name,
                datastore_filename='zinc10_sample.csv',
                filename='./assets/zinc10.csv',
                description='Sample test csv file')
    
    featurized_file_address = featurize(profile_name,
              project_name,
              dataset_address,
              featurizer='ecfp',
              output='test_featurized',
              dataset_column='smiles',
              label_column='logp',
              feat_kwargs = {'size': 1024},
              )
    
    assert featurized_file_address == f'deepchem://{profile_name}/{project_name}/test_featurized'
