import os
from pathlib import Path
import re
import tempfile

import pandas as pd  # type: ignore
from pandas.testing import assert_frame_equal  # type: ignore
import pytest

from deepchem_server.core import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.primitives.fep.rbfe.collate_rbfe_results import (
    calculate_dg_values,
    collate_rbfe_results,
    get_ligands_from_results,
    get_result_dataframe,
    process_input_files,
)


ASSETS_DIR = Path(__file__).resolve().parent.parent.parent / "assets"


@pytest.fixture()
def valid_uploaded_files_addresses():
    VALID_RESULTS_DIR = ASSETS_DIR / "rbfe_results" / "valid_results"

    uploaded_files_deepchem_addresses = []

    for item in VALID_RESULTS_DIR.iterdir():
        if item.is_file():
            item_datacard = DataCard(address='', file_type='json', data_type='text/plain')
            address = config.get_datastore().upload_data(item.name, item.as_posix(), item_datacard)
            uploaded_files_deepchem_addresses.append(address)
    return uploaded_files_deepchem_addresses


def test_process_input_files_valid_input(valid_uploaded_files_addresses):

    simulation_results = process_input_files(valid_uploaded_files_addresses)

    assert len(simulation_results) == 4

    for result in simulation_results:
        assert result['componentA_name'] is not None
        assert result['componentB_name'] is not None
        assert result['edge_ddG'] is not None
        assert result['edge_ddG_uncertainty'] is not None


def test_process_input_files_with_invalid_file(caplog):
    VALID_AND_INVALID_RESULTS = ASSETS_DIR / "rbfe_results" / "valid_and_invalid_results"

    uploaded_files_deepchem_addresses = []

    for item in VALID_AND_INVALID_RESULTS.iterdir():
        if item.is_file():
            item_datacard = DataCard(address='', file_type='json', data_type='text/plain')
            address = config.get_datastore().upload_data(item.name, item.as_posix(), item_datacard)
            uploaded_files_deepchem_addresses.append(address)

    simulation_results = process_input_files(uploaded_files_deepchem_addresses)

    assert len(simulation_results) == 3

    assert re.match(r".*Could not process(.+?){}.*$".format(re.escape("rbfe-lig_ejm_48-lig_jmc_28.json")), caplog.text)


def test_get_ligands_from_results(valid_uploaded_files_addresses):
    simulation_results = process_input_files(valid_uploaded_files_addresses)

    ligands = get_ligands_from_results(simulation_results)

    assert type(ligands) is list
    assert len(ligands) == 5
    assert sorted(ligands) == ['lig_ejm_31', 'lig_ejm_46', 'lig_ejm_48', 'lig_jmc_23', 'lig_jmc_28']


def test_calculate_dg_values(valid_uploaded_files_addresses):
    simulation_results = process_input_files(valid_uploaded_files_addresses)

    ligand_dg_values = calculate_dg_values(simulation_results, "lig_ejm_48", "-9.00 kcal/mol", "0.00 kcal/mol")

    rounded_dg_values = [round(dg_value.value.magnitude, 2) for dg_value in ligand_dg_values.values()]
    assert sorted(rounded_dg_values) == sorted([-10.72, -9.93, -10.39, -9.0, -10.78])

    rounded_dg_uncertainties = [round(dg_value.uncertainty.magnitude, 3) for dg_value in ligand_dg_values.values()]
    assert sorted(rounded_dg_uncertainties) == [0.0, 0.105, 0.16, 0.165, 0.421]


def test_get_result_dataframe(caplog, valid_uploaded_files_addresses):
    simulation_results = process_input_files(valid_uploaded_files_addresses)

    ligand_dg_values = calculate_dg_values(simulation_results, "lig_ejm_48", "-9.00 kcal/mol", "0.00 kcal/mol")

    result_dataframe = get_result_dataframe(ligand_dg_values)

    expected_dataframe = pd.read_csv(ASSETS_DIR / "rbfe_results" / "expected_dataframe.csv")

    result_dataframe = result_dataframe.astype('string')
    expected_dataframe = expected_dataframe.astype('string')

    assert_frame_equal(result_dataframe, expected_dataframe)


def test_collate_rbfe_results(valid_uploaded_files_addresses):
    output_file_name = "pytest_output.csv"
    output_datastore_address = collate_rbfe_results(valid_uploaded_files_addresses, 'lig_ejm_48', '-9.00 kcal/mol',
                                                    '0.00 kcal/mol', output_file_name)

    tempdir = tempfile.TemporaryDirectory()
    output_local_address = os.path.join(tempdir.name, 'downloaded_csv.csv')
    config.get_datastore().download_object(output_datastore_address, output_local_address)

    # Assert that the csv file is present in the temporary directory
    assert os.path.exists(output_local_address)

    result_dataframe = pd.read_csv(output_local_address)
    result_dataframe = result_dataframe.astype('string')

    expected_dataframe = pd.read_csv(ASSETS_DIR / "rbfe_results" / "expected_dataframe.csv")

    expected_dataframe = expected_dataframe.astype('string')

    assert_frame_equal(result_dataframe, expected_dataframe)
