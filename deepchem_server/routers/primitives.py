import json
from typing import Dict, Optional

from fastapi import APIRouter, HTTPException

from deepchem_server.core import model_mappings
from deepchem_server.core.feat import featurizer_map
from deepchem_server.utils import parse_boolean_none_values_from_kwargs, run_job

router = APIRouter(
    prefix="/primitive",
    tags=["primitive"],
)


@router.post("/featurize")
async def featurize(
    profile_name: str,
    project_name: str,
    dataset_address: str,
    featurizer: str,
    output: str,
    dataset_column: str,
    feat_kwargs: Optional[Dict] = None,
    label_column: Optional[str] = None,
) -> dict:
    """
    Submits a featurization job

    Parameters
    ----------
    profile_name: str
        Name of the Profile where the job is run
    project_name: str
        Name of the Project where the job is run
    dataset_address: str
        datastore address of dataset to featurize
    featurizer: str
        featurizer to use
    output: str
        name of the featurized dataset
    dataset_column: Optional[str]
        Column containing the input for featurizer
    feat_kwargs: Optional[Dict]
        Keyword arguments to pass to featurizer on initialization
    label_column: Optional[str]
        The target column in case this dataset is going to be used for training purposes

    Raises
    ------
    HTTPException
        If the featurizer is invalid or if the featurization fails.

    Returns
    -------
    dict
        A dictionary containing the address of the featurized dataset.
    """
    if not feat_kwargs:
        feat_kwargs = {'feat_kwargs': {}}

    if featurizer not in featurizer_map.keys():
        message = f"Invalid featurizer: {featurizer}. Use one of {list(featurizer_map.keys())}."
        raise HTTPException(status_code=404, detail=message)

    if isinstance(feat_kwargs['feat_kwargs'], str):
        feat_kwargs['feat_kwargs'] = json.loads(feat_kwargs['feat_kwargs'])

    for key, value in feat_kwargs['feat_kwargs'].items():
        if isinstance(value, str):
            if value.lower() == "true":
                feat_kwargs['feat_kwargs'][key] = True
            elif value.lower() == "false":
                feat_kwargs['feat_kwargs'][key] = False
            elif value.lower() == "none":
                feat_kwargs['feat_kwargs'][key] = None

    program: Dict = {
        'program_name': 'featurize',
        'dataset_address': dataset_address,
        'featurizer': featurizer,
        'output': output,
        'dataset_column': dataset_column,
        'feat_kwargs': feat_kwargs['feat_kwargs'],
        'label_column': label_column,
    }

    try:
        result = run_job(profile_name=profile_name,
                         project_name=project_name,
                         program=program)
    except Exception as e:
        raise HTTPException(status_code=500,
                            detail=f"Featurization failed: {str(e)}")

    return {"featurized_file_address": str(result)}


@router.post("/train")
async def train(
    profile_name: str,
    project_name: str,
    dataset_address: str,
    model_type: str,
    model_name: str,
    init_kwargs: Optional[Dict] = None,
    train_kwargs: Optional[Dict] = None,
) -> dict:
    """
    Submits a training job

    Parameters
    ----------
    profile_name: str
        Name of the Profile where the job is run
    project_name: str
        Name of the Project where the job is run
    dataset_address: str
        datastore address of dataset to train on
    model_type: str
        type of model to train
    model_name: str
        name of the trained model
    init_kwargs: Optional[Dict]
        Keyword arguments to pass to model on initialization
    train_kwargs: Optional[Dict]
        Keyword arguments to pass to model on training

    Raises
    ------
    HTTPException
        If the model type is invalid or if the training fails.

    Returns
    -------
    dict
        A dictionary containing the address of the trained model.
    """
    if not init_kwargs:
        init_kwargs = {}
    if not train_kwargs:
        train_kwargs = {}

    if model_type not in model_mappings.model_address_map.keys():
        message = f"Invalid model type: {model_type}. Use one of {list(model_mappings.model_address_map.keys())}."
        raise HTTPException(status_code=404, detail=message)

    if isinstance(init_kwargs, str):
        init_kwargs = json.loads(init_kwargs)

    if isinstance(train_kwargs, str):
        train_kwargs = json.loads(train_kwargs)

    init_kwargs = parse_boolean_none_values_from_kwargs(init_kwargs)
    train_kwargs = parse_boolean_none_values_from_kwargs(train_kwargs)

    program: Dict = {
        "program_name": "train",
        "dataset_address": dataset_address,
        "model_type": model_type,
        "model_name": model_name,
        "init_kwargs": init_kwargs,
        "train_kwargs": train_kwargs,
    }

    try:
        result = run_job(profile_name=profile_name,
                         project_name=project_name,
                         program=program)
    except Exception as e:
        raise HTTPException(status_code=500,
                            detail=f"Training failed: {str(e)}")

    return {"trained_model_address": str(result)}


@router.post("/evaluate")
async def evaluate(
    profile_name: str,
    project_name: str,
    dataset_addresses: List[str],
    model_address: str,
    metrics: List[str],
    output_key: str,
    is_metric_plots: bool = False,
) -> Union[Dict, JSONResponse]:
    """
    Submits an evaluation job

    Parameters
    ----------
    profile_name: str
        Name of the Profile where the job is run
    project_name: str
        Name of the Project where the job is run
    dataset_addresses: List[str]
        List of dataset addresses to evaluate the model on
    model_address: str
        datastore address of the trained model
    metrics: List[str]
        List of metrics to evaluate the model with
    output_key: str
        Name of the evaluation output
    is_metric_plots: bool
        Whether plot based metric is used or not
    """

    program: Dict = {
        "program_name": "evaluate",
        "dataset_addresses": dataset_addresses,
        "model_address": model_address,
        "metrics": metrics,
        "output_key": output_key,
        "is_metric_plots": is_metric_plots,
    }

    result = run_job(profile_name=profile_name,
                     project_name=project_name,
                     program=program)

    if isinstance(result, Exception):
        return JSONResponse(
            status_code=500,
            content={"message": f"Evaluation failed: {str(result)}"})

    return {"evaluation_result_address": str(result)}
