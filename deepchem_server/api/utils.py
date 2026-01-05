"""
API utilities for request parsing and data handling.

This module contains utility functions used by API routers for:
- Uploading data to the datastore
- Parsing request parameters
"""
import ast
import json
import logging
from typing import Dict

from deepchem_server.core.common import config


logger = logging.getLogger(__name__)


def upload_data(profile_name, project_name, datastore_filename, contents, data_card, backend='local'):
    """
    Upload data to datastore.
    
    Uses DatastoreClient to communicate with the datastore service.

    Parameters
    ----------
    profile_name: str
        Name of the Profile where the job is run
    project_name: str
        Name of the Project where the job is run
    datastore_filename: str
        The file name in which data is to be written in DataStore
    contents: bytes
        The file content to upload
    data_card: DataCard
        data card for the file
    backend: str
        Backend to be used (ignored, kept for compatibility)
        
    Returns
    -------
    str
        The deepchem address of the uploaded data
    """
    client = config.get_datastore_client()
    if client is None:
        raise RuntimeError("DATASTORE_URL environment variable must be set. "
                           "All operations require a running datastore service.")

    address = f"deepchem://{profile_name}/{project_name}/{datastore_filename}"
    data_card.address = address
    card_dict = json.loads(data_card.to_json())

    return client.upload_data(address, contents, card_dict, kind="data")


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


def parse_dict_with_datatypes(dict_to_parse):
    """
    Parses a dictionary with string values and returns a new dictionary with
    the same keys but with the values converted to their corresponding datatypes.

    Parameters
    ----------
    dict_to_parse: dict
        The dictionary to parse.

    Returns
    -------
    parsed_dict: dict
        A new dictionary with the same keys as the input dictionary but with
        the values converted to their corresponding datatypes.
    """
    parsed_dict = {}
    if dict_to_parse is None:
        return parsed_dict

    for key, value in dict_to_parse.items():
        # Check if the value is a non-empty string
        parsed_value = value
        if isinstance(value, str):
            # Evaluate the string using ast.literal_eval(). Throws ValueError and SyntaxError on failure
            if value.lower() == "true":
                parsed_value = True
            elif value.lower() == "false":
                parsed_value = False
            elif value.lower() == "none":
                parsed_value = None
            elif value and len(value) < 10000:
                parsed_value = ast.literal_eval(value)
        parsed_dict[key] = parsed_value
    return parsed_dict
