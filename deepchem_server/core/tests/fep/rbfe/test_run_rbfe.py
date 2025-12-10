import os
import tempfile  # noqa
from unittest.mock import patch

from openfe import SolventComponent  # type: ignore
from openff.units import unit  # type: ignore
import pint
import pytest

from deepchem_server.core import config
from deepchem_server.core.primitives.fep.rbfe.data_domain_classes.EdgeSimulationResult import EdgeSimulationResult
from deepchem_server.core.primitives.fep.rbfe.run_rbfe import run_rbfe
from deepchem_server.core.primitives.fep.rbfe.utils.constants import NetworkPlanningConstants


mocked_result = EdgeSimulationResult(
    componentA_name='benzene',
    componentB_name="toluene",
    complex_dG=pint.Quantity('1 kcal/mole'),  # type: ignore
    is_dry_run=False,
    dry_run_status=None,
    complex_dG_uncertainty=pint.Quantity('1 kcal/mole'),  # type: ignore
    solvent_dG=pint.Quantity('1 kcal/mole'),  # type: ignore
    solvent_dG_uncertainty=pint.Quantity('1 kcal/mole'),  # type: ignore
    edge_ddG=pint.Quantity('1 kcal/mole'),  # type: ignore
    edge_ddG_uncertainty=pint.Quantity('1 kcal/mole'),  # type: ignore
)


@pytest.fixture()
def run_rbfe_kwargs(ligands_large_datastore_address, protein_datastore_address):

    solvent = SolventComponent(positive_ion="Na",
                               negative_ion="Cl",
                               ion_concentration=0.15 * unit.molar,
                               neutralize=True)

    solvent = '{"positive_ion": "Na", "negative_ion": "Cl", \
        "ion_concentration": "0.15 mol/l", "neutralize": True}'

    # In openFE v1.1.0: n_steps = time_per_iteration/timestep
    # In openFE v1.1.0: alchemical_sampler_settings.n_repeats is moved to rbfe_settings.protocol_repeats
    overridden_rbfe_settings = '{"protocol_repeats": "1", \
        "simulation_settings" : {"equilibration_length" : "0.1 ps", \
            "production_length" : "0.1 ps", "time_per_iteration" : "0.01 ps"},\
                  "integrator_settings": {"timestep": "0.01 ps"}}'

    return {
        "ligands_sdf_address": ligands_large_datastore_address,
        "cleaned_protein_pdb_address": protein_datastore_address,
        "network_type": NetworkPlanningConstants.PerturbationNetworkType.MINIMAL_SPANNING.name,
        "scorer_type": "LOMAP",
        "solvent_json": solvent,
        "overridden_rbfe_settings": overridden_rbfe_settings,
        "radial_network_central_ligand": "benzene"
    }


@pytest.fixture()
def run_rbfe_nested_kwargs(ligands_nested_large_datastore_address, protein_nested_datastore_address):

    solvent = SolventComponent(positive_ion="Na",
                               negative_ion="Cl",
                               ion_concentration=0.15 * unit.molar,
                               neutralize=True)

    solvent = '{"positive_ion": "Na", "negative_ion": "Cl", \
        "ion_concentration": "0.15 mol/l", "neutralize": True}'

    # In openFE v1.1.0: n_steps = time_per_iteration/timestep
    # In openFE v1.1.0: alchemical_sampler_settings.n_repeats is moved to rbfe_settings.protocol_repeats
    overridden_rbfe_settings = '{"simulation_settings" : \
        {"equilibration_length" : "0.1 ps", "production_length" : "0.1 ps", \
            "time_per_iteration" : "0.01 ps"}, "integrator_settings": \
                {"timestep": "0.01 ps"}}'

    return {
        "ligands_sdf_address": ligands_nested_large_datastore_address,
        "cleaned_protein_pdb_address": protein_nested_datastore_address,
        "network_type": NetworkPlanningConstants.PerturbationNetworkType.MINIMAL_SPANNING.name,
        "scorer_type": "LOMAP",
        "solvent_json": solvent,
        "overridden_rbfe_settings": overridden_rbfe_settings,
    }


@patch('deepchem_server.core.primitives.fep.rbfe.run_rbfe.system_setup.dry_run_edge', return_value=mocked_result)
def test_run_rbfe_dry_run(patched_run_edge, caplog, protein_datastore_address, ligands_large_datastore_address,
                          run_rbfe_kwargs):
    """Tests RBFE dry run
    """

    # Without the AWS Batch environment variable
    try:
        del os.environ["AWS_BATCH_JOB_ARRAY_INDEX"]
    except KeyError:
        pass

    result_datastore_addresses = run_rbfe(**run_rbfe_kwargs, dry_run=True)

    assert result_datastore_addresses is not None

    # Assert that the log message for AWS Batch confirms unavailability of AWS Batch Environment.
    assert 'Could not find AWS_BATCH_JOB_ARRAY_INDEX. Running all' in caplog.text
    assert len(result_datastore_addresses) == 6

    tempdir = tempfile.TemporaryDirectory()

    for index, result_datastore_address in enumerate(result_datastore_addresses):

        assert result_datastore_address is not None

        result_destination_address = os.path.join(tempdir.name, f'result-{index}.json')
        config.get_datastore().download_object(result_datastore_address, result_destination_address)

        # Assert that the file result.json is present in the temporary directory
        assert os.path.exists(result_destination_address)

    # With the AWS Batch environment variable
    os.environ['AWS_BATCH_JOB_ARRAY_INDEX'] = '3'
    result_datastore_addresses = run_rbfe(**run_rbfe_kwargs, dry_run=True)

    # Assert that the log message for AWS Batch confirms availability of AWS Batch Environment.
    assert "Found AWS_BATCH_JOB_ARRAY_INDEX in the environment. Running edge 3" in caplog.text


@patch('deepchem_server.core.primitives.fep.rbfe.run_rbfe.system_setup.dry_run_edge', return_value=mocked_result)
def test_run_rbfe_dry_run_nested_full_address(patched_run_edge, caplog, protein_nested_datastore_address,
                                              ligands_nested_large_datastore_address, run_rbfe_nested_kwargs):
    """Tests RBFE dry run
    """

    # Without the AWS Batch environment variable
    try:
        del os.environ["AWS_BATCH_JOB_ARRAY_INDEX"]
    except KeyError:
        ...

    result_datastore_addresses = run_rbfe(**run_rbfe_nested_kwargs, dry_run=True)

    assert result_datastore_addresses is not None
    assert 'deepchem://test/user/test rbfe/rbfe-benzene-phenol.json' in result_datastore_addresses


@patch('deepchem_server.core.primitives.fep.rbfe.run_rbfe.system_setup.run_edge', return_value=mocked_result)
def test_run_rbfe(patched_run_edge, caplog, protein_datastore_address, ligands_large_datastore_address,
                  run_rbfe_kwargs):
    """Tests RBFE dry run
    """

    # Without the AWS Batch environment variable
    try:
        del os.environ["AWS_BATCH_JOB_ARRAY_INDEX"]
    except KeyError:
        pass

    result_datastore_addresses = run_rbfe(**run_rbfe_kwargs, dry_run=False)

    assert result_datastore_addresses is not None

    # Assert that the log message for AWS Batch confirms unavailability of AWS Batch Environment.
    assert 'Could not find AWS_BATCH_JOB_ARRAY_INDEX. Running all' in caplog.text
    assert len(result_datastore_addresses) == 6

    tempdir = tempfile.TemporaryDirectory()

    for index, result_datastore_address in enumerate(result_datastore_addresses):

        assert result_datastore_address is not None

        result_destination_address = os.path.join(tempdir.name, f'result-{index}.json')
        config.get_datastore().download_object(result_datastore_address, result_destination_address)

        # Assert that the file result.json is present in the temporary directory
        assert os.path.exists(result_destination_address)

    # With the AWS Batch environment variable
    os.environ['AWS_BATCH_JOB_ARRAY_INDEX'] = '3'
    result_datastore_addresses = run_rbfe(**run_rbfe_kwargs, dry_run=False)

    # Assert that the log message for AWS Batch confirms availability of AWS Batch Environment.
    assert "Found AWS_BATCH_JOB_ARRAY_INDEX in the environment. Running edge 3" in caplog.text


@patch('deepchem_server.core.primitives.fep.rbfe.run_rbfe.system_setup.run_edge', return_value=mocked_result)
def test_run_rbfe_single_dir_output_key(patched_run_edge,
                                        caplog,
                                        protein_datastore_address,
                                        ligands_large_datastore_address,
                                        run_rbfe_kwargs: dict,
                                        output_key: str = "output_key"):
    """Tests RBFE dry run, with specific output key, where the output key is a single directory.

    Args:
        run_rbfe_kwargs (dict): Args for running run_rbfe() function.
        output_key (str, optional): Folder where the data is going to get stored. Defaults to "output_key".
    """

    # Without the AWS Batch environment variable
    try:
        del os.environ["AWS_BATCH_JOB_ARRAY_INDEX"]
    except KeyError:
        pass

    result_datastore_addresses = run_rbfe(**run_rbfe_kwargs, dry_run=False, output_key=output_key)

    assert result_datastore_addresses is not None

    # Assert that the log message for AWS Batch confirms unavailability of AWS Batch Environment.
    assert 'Could not find AWS_BATCH_JOB_ARRAY_INDEX. Running all' in caplog.text
    assert len(result_datastore_addresses) == 6

    # Check presnece of output_key
    for address in result_datastore_addresses:
        assert output_key in address, "output_key not in address"

    # With the AWS Batch environment variable
    os.environ['AWS_BATCH_JOB_ARRAY_INDEX'] = '3'
    result_datastore_addresses = run_rbfe(**run_rbfe_kwargs, dry_run=False, output_key=output_key)

    # Assert that the log message for AWS Batch confirms availability of AWS Batch Environment.
    assert "Found AWS_BATCH_JOB_ARRAY_INDEX in the environment. Running edge 3" in caplog.text
    assert output_key in result_datastore_addresses[0], \
        "output_key not in address"


@patch('deepchem_server.core.primitives.fep.rbfe.run_rbfe.system_setup.run_edge', return_value=mocked_result)
def test_run_rbfe_nested_dir_output_key(patched_run_edge,
                                        caplog,
                                        protein_datastore_address,
                                        ligands_large_datastore_address,
                                        run_rbfe_kwargs: dict,
                                        output_key: str = "output_key/key_2/key_3"):
    """Tests RBFE dry run, with specific output key, where the output key is a nested directory.

    Args:
        run_rbfe_kwargs (dict): Args for running run_rbfe() function.
        output_key (str, optional): Folder(here nested folder) where the data is going to get stored. Defaults to "output_key/key_2/key_3".
    """

    # Without the AWS Batch environment variable
    try:
        del os.environ["AWS_BATCH_JOB_ARRAY_INDEX"]
    except KeyError:
        pass

    result_datastore_addresses = run_rbfe(**run_rbfe_kwargs, dry_run=False, output_key=output_key)

    assert result_datastore_addresses is not None

    # Assert that the log message for AWS Batch confirms unavailability of AWS Batch Environment.
    assert 'Could not find AWS_BATCH_JOB_ARRAY_INDEX. Running all' in caplog.text
    assert len(result_datastore_addresses) == 6

    # Check presnece of output_key
    for address in result_datastore_addresses:
        assert output_key in address, "output_key not in address"

    # With the AWS Batch environment variable
    os.environ['AWS_BATCH_JOB_ARRAY_INDEX'] = '3'
    result_datastore_addresses = run_rbfe(**run_rbfe_kwargs, dry_run=False, output_key=output_key)

    # Assert that the log message for AWS Batch confirms availability of AWS Batch Environment.
    assert "Found AWS_BATCH_JOB_ARRAY_INDEX in the environment. Running edge 3" in caplog.text
    assert output_key in result_datastore_addresses[0], \
        "output_key not in address"
