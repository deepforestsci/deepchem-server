from deepchem_server.core.primitives.fep.rbfe import system_setup
import pytest  # type: ignore
from openfe import (  # type: ignore
    ProteinComponent, SmallMoleculeComponent, SolventComponent,
)
from deepchem_server.core.primitives.fep.rbfe.data_domain_classes.RunnableEdge import RunnableEdge
from openff.units import unit  # type: ignore
from gufe.ligandnetwork import LigandNetwork  # noqa
from deepchem_server.core.primitives.fep.rbfe.data_domain_classes.EdgeSimulationResult import EdgeSimulationResult
from deepchem_server.core.primitives.fep.rbfe.utils.constants import EMPTY_STRING, SUCCESS, TOLUENE, BENZENE, NetworkPlanningConstants  # noqa


def test_load_ligands(ligands_datastore_address):
    ligand_mols = system_setup.load_ligands(ligands_datastore_address)

    assert isinstance(ligand_mols, list)

    assert len(ligand_mols) == 2

    for ligand in ligand_mols:
        assert isinstance(ligand, SmallMoleculeComponent)


def test_get_reference_ligand(ligands_datastore_address):
    ligand_mols = system_setup.load_ligands(ligands_datastore_address)

    reference_ligand_index = system_setup.get_reference_ligand(ligand_mols,)

    assert ligand_mols[reference_ligand_index].name == BENZENE


def test_get_perturbation_network_invalid_input(ligands_datastore_address):

    ligand_mols = system_setup.load_ligands(ligands_datastore_address)

    with pytest.raises(AssertionError):
        system_setup.get_perturbation_network(None, NetworkPlanningConstants.PerturbationNetworkType.RADIAL, None)

    with pytest.raises(AssertionError):
        system_setup.get_perturbation_network(ligand_mols, None, None)


def test_get_perturbation_network_no_central_ligand_for_radial_network(ligands_datastore_address):

    ligand_mols = system_setup.load_ligands(ligands_datastore_address)

    network = system_setup.get_perturbation_network(ligand_mols,
                                                    NetworkPlanningConstants.PerturbationNetworkType.RADIAL)

    # Assert that central ligand is Benzene
    assert list(network.edges)[0].componentA.name == BENZENE


def test_get_perturbation_network_no_scorer_for_maximal_and_minimal(ligands_datastore_address):

    ligand_mols = system_setup.load_ligands(ligands_datastore_address)

    with pytest.raises(AssertionError):
        system_setup.get_perturbation_network(ligand_mols,
                                              NetworkPlanningConstants.PerturbationNetworkType.MINIMAL_SPANNING)

    with pytest.raises(AssertionError):
        system_setup.get_perturbation_network(ligand_mols, NetworkPlanningConstants.PerturbationNetworkType.MAXIMAL)


def test_get_perturbation_network_maximal(ligands_large_datastore_address):
    ligand_mols = system_setup.load_ligands(ligands_large_datastore_address)

    network = system_setup.get_perturbation_network(ligand_mols,
                                                    NetworkPlanningConstants.PerturbationNetworkType.MAXIMAL,
                                                    NetworkPlanningConstants.ScorerType.LOMAP)

    assert len(list(network.edges)) == 27


def test_get_perturbation_network_minimal_spanning(ligands_large_datastore_address):
    ligand_mols = system_setup.load_ligands(ligands_large_datastore_address)

    network = system_setup.get_perturbation_network(ligand_mols,
                                                    NetworkPlanningConstants.PerturbationNetworkType.MINIMAL_SPANNING,
                                                    NetworkPlanningConstants.ScorerType.LOMAP)

    assert len(network.edges) == 6


def test_load_cleaned_protein(protein_datastore_address):

    protein = system_setup.load_cleaned_protein(protein_datastore_address)

    assert isinstance(protein, ProteinComponent)


def test_default_settings():

    # Load the default RBFE settings
    default_settings = system_setup.get_default_RBFE_simulation_settings()

    assert default_settings.simulation_settings.equilibration_length == 2 * unit.nanosecond
    assert default_settings.simulation_settings.production_length == 10 * unit.nanosecond
    assert default_settings.simulation_settings.real_time_analysis_interval is None


def test_setup_transformations(ligands_datastore_address, protein_datastore_address):
    ligand_mols = system_setup.load_ligands(ligands_datastore_address)

    # Define the transformation network.
    network = system_setup.get_perturbation_network(ligand_mols,
                                                    NetworkPlanningConstants.PerturbationNetworkType('RADIAL'),
                                                    NetworkPlanningConstants.ScorerType('LOMAP'),
                                                    central_ligand="benzene")

    solvent = SolventComponent(positive_ion='Na', negative_ion='Cl', ion_concentration=0.15 * unit.molar)
    cleaned_protein = system_setup.load_cleaned_protein(protein_datastore_address)
    rbfe_settings = system_setup.get_default_RBFE_simulation_settings()

    runnable_edges, rbfe_transform = system_setup.setup_transformations(network, solvent, cleaned_protein,
                                                                        rbfe_settings)

    assert isinstance(runnable_edges, list)

    for runnable_edge in runnable_edges:
        assert isinstance(runnable_edge, RunnableEdge)

    assert rbfe_transform is not None


def test_dry_run_edge(ligands_datastore_address, protein_datastore_address):
    ligand_mols = system_setup.load_ligands(ligands_datastore_address)

    network = system_setup.get_perturbation_network(
        ligand_mols,
        NetworkPlanningConstants.PerturbationNetworkType('MINIMAL_SPANNING'),
        NetworkPlanningConstants.ScorerType('LOMAP'),
        central_ligand='benzene')

    solvent = SolventComponent(positive_ion='Na', negative_ion='Cl', ion_concentration=0.15 * unit.molar)
    cleaned_protein = system_setup.load_cleaned_protein(protein_datastore_address)

    rbfe_settings = system_setup.get_default_RBFE_simulation_settings()

    # Tweak settings for lesser run times.
    rbfe_settings.simulation_settings.equilibration_length = 0.1 * unit.picosecond
    rbfe_settings.simulation_settings.production_length = 0.1 * unit.picosecond
    rbfe_settings.simulation_settings.time_per_iteration = 0.01 * unit.picosecond
    rbfe_settings.integrator_settings.timestep = 0.01 * unit.picosecond
    rbfe_settings.protocol_repeats = 1
    rbfe_settings.output_settings.checkpoint_interval = 5.0 * unit.picosecond
    rbfe_settings.engine_settings.compute_platform = 'CPU'

    runnable_edges, _ = system_setup.setup_transformations(network, solvent, cleaned_protein, rbfe_settings)

    edgeSimulationResult = system_setup.dry_run_edge(runnable_edges[0])
    assert isinstance(edgeSimulationResult, EdgeSimulationResult)

    # Expected attributes should be present
    assert edgeSimulationResult.is_dry_run is True
    assert edgeSimulationResult.dry_run_status == SUCCESS
    assert edgeSimulationResult.componentA_name is not None
    assert edgeSimulationResult.componentB_name is not None

    # The following attributes should be empty in a Dry Run
    assert edgeSimulationResult.complex_dG == EMPTY_STRING
    assert edgeSimulationResult.complex_dG_uncertainty == EMPTY_STRING
    assert edgeSimulationResult.solvent_dG == EMPTY_STRING
    assert edgeSimulationResult.solvent_dG_uncertainty == EMPTY_STRING
    assert edgeSimulationResult.edge_ddG == EMPTY_STRING
    assert edgeSimulationResult.edge_ddG_uncertainty == EMPTY_STRING
