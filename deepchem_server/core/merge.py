"""This file contains utilities to merge datasets."""

import ast
import csv
import os
import pathlib
import tempfile
from typing import List

from deepchem_server.core import config
from deepchem_server.core.address import DeepchemAddress
from deepchem_server.core.cards import DataCard


def _merge_csv(dataset_addresses: List[str], output_key: str):
    """Merge multiple csv files into a single csv file.
    Parameters
    ----------
    dataset_addresses : List[str]
        List of dataset addresses to merge.
    output_key : str
        The output key for the merged dataset.

    Returns
    -------
    str
        The address of the merged dataset.
    """
    datastore = config.get_datastore()
    input_file_paths = []
    tempdir_input_files = tempfile.TemporaryDirectory()
    for i, address in enumerate(dataset_addresses):
        file_download_path = pathlib.Path(tempdir_input_files.name) / str(i)
        # Convert path to string
        input_file_paths.append(str(file_download_path))
        if datastore is None:
            raise ValueError("Datastore not set")
        datastore.download_object(address, str(file_download_path))

    input_file = input_file_paths[0]
    with open(input_file, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        fieldnames = reader.fieldnames
    del input_file

    # Ensure fieldnames is not None
    if fieldnames is None:
        raise ValueError("Fieldnames cannot be None")

    tempdir = tempfile.TemporaryDirectory()
    output_key = DeepchemAddress.get_key(output_key)
    temp_output_path = os.path.join(tempdir.name, output_key)
    temp_output_dir = os.path.dirname(temp_output_path)
    if not os.path.exists(temp_output_dir):
        os.makedirs(temp_output_dir)
    with open(temp_output_path, "w") as outputfile:
        writer = csv.DictWriter(outputfile, fieldnames=fieldnames)
        writer.writeheader()

        for file in input_file_paths:
            with open(file, "r") as inputfile:
                reader = csv.DictReader(inputfile)
                for row in reader:
                    writer.writerow(row)

    card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    if datastore is None:
        raise ValueError("Datastore not set")
    address = datastore.upload_data(output_key, temp_output_path, card)
    return address


def merge(dataset_addresses: List[str], output_key: str):
    """Merge multiple datasets into a single dataset.
    Parameters
    ----------
    dataset_addresses : List[str]
        List of dataset addresses to merge.
    output_key : str
        The output key for the merged dataset.

    Returns
    -------
    str
        The address of the merged dataset.
    """
    if isinstance(dataset_addresses, str):
        dataset_addresses = ast.literal_eval(dataset_addresses)
    if output_key.endswith(".csv"):
        return _merge_csv(dataset_addresses, output_key)
    else:
        raise TypeError("merge operation is only available for csv file types")
