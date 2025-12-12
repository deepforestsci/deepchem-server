from deepchem_server.core.primitives.fep.rbfe.utils.rbfe_utils import RBFESettingsUtils
from openff.units import unit  # type: ignore
import pytest  # type: ignore
from openfe import SolventComponent  # type: ignore
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol  # type: ignore


def test_rbfe_utils_loads_NoneType_input():
    with pytest.raises(TypeError):
        RBFESettingsUtils.loads(None)


def test_rbfe_utils_loads_invalid_input():

    # Invalid JSON string. The correct key is "simulation_settings", and not "simulation_setting".
    invalid_rbfe_settings = '{"protocol_repeats": 1, "simulation_setting" : {"equilibration_length" : "2 ns", "production_length" : "2 ns"}}'
    with pytest.raises(AttributeError):
        RBFESettingsUtils.loads(invalid_rbfe_settings)


def test_rbfe_utils_loads_valid_input():

    # A valid rbfe_settings JSON string with two overridden parameters.
    overridden_rbfe_settings = '{"simulation_settings" : {"equilibration_length" : "2 ns", "production_length" : "20 ns"}}'

    # Loading the JSON
    rbfe_settings = RBFESettingsUtils.loads(overridden_rbfe_settings)

    assert rbfe_settings.simulation_settings.equilibration_length == 2 * unit.nanosecond
    assert rbfe_settings.simulation_settings.production_length == 20 * unit.nanosecond


def test_rbfe_utils_dumps_NoneType_input():
    with pytest.raises(TypeError):
        RBFESettingsUtils.dumps(None)


def test_rbfe_utils_dumps_invalid_input():
    # Testing with a SolventComponent object
    incorrect_type_object = SolventComponent()

    with pytest.raises(TypeError):
        RBFESettingsUtils.dumps(incorrect_type_object)


def test_rbfe_utils_dumps_valid_input():

    # Initializing a default RBFESettings object
    rbfe_settings = RelativeHybridTopologyProtocol.default_settings()

    # Serializing the RBFESettings object
    rbfe_settings_json = RBFESettingsUtils.dumps(rbfe_settings)

    assert rbfe_settings_json == '{"forcefield_settings": {"constraints": "hbonds", "rigid_water": true, "hydrogen_mass": 3.0, "forcefields": ["amber/ff14SB.xml", "amber/tip3p_standard.xml", "amber/tip3p_HFE_multivalent.xml", "amber/phosaa10.xml"], "small_molecule_forcefield": "openff-2.1.1", "nonbonded_cutoff": "1.0 nanometer", "nonbonded_method": "PME"}, "thermo_settings": {"temperature": "298.15 kelvin", "pressure": "0.9869232667160129 standard_atmosphere", "ph": null, "redox_potential": null}, "protocol_repeats": 3, "solvation_settings": {"solvent_model": "tip3p", "solvent_padding": "1.2 nanometer", "box_shape": "cube", "number_of_solvent_molecules": null, "box_vectors": null, "box_size": null}, "partial_charge_settings": {"partial_charge_method": "am1bcc", "off_toolkit_backend": "ambertools", "number_of_conformers": null, "nagl_model": null}, "lambda_settings": {"lambda_functions": "default", "lambda_windows": 11}, "alchemical_settings": {"softcore_LJ": "gapsys", "explicit_charge_correction_cutoff": "0.8 nanometer", "endstate_dispersion_correction": false, "use_dispersion_correction": false, "softcore_alpha": 0.85, "turn_off_core_unique_exceptions": false, "explicit_charge_correction": false}, "simulation_settings": {"equilibration_length": "1.0 nanosecond", "production_length": "5.0 nanosecond", "minimization_steps": 5000, "time_per_iteration": "1 picosecond", "real_time_analysis_interval": "250 picosecond", "early_termination_target_error": "0.0 kilocalorie_per_mole", "real_time_analysis_minimum_time": "500 picosecond", "sampler_method": "repex", "sams_flatness_criteria": "logZ-flatness", "sams_gamma0": 1.0, "n_replicas": 11}, "engine_settings": {"compute_platform": null}, "integrator_settings": {"timestep": "4 femtosecond", "langevin_collision_rate": "1.0 / picosecond", "barostat_frequency": "25 timestep", "remove_com": false, "reassign_velocities": false, "n_restart_attempts": 20, "constraint_tolerance": 1e-06}, "output_settings": {"checkpoint_interval": "250 picosecond", "forcefield_cache": "db.json", "output_indices": "not water", "checkpoint_storage_filename": "checkpoint.chk", "output_filename": "simulation.nc", "output_structure": "hybrid_system.pdb"}}'  # noqa: E501
