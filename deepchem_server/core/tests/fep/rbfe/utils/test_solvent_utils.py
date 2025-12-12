from deepchem_server.core.primitives.fep.rbfe.utils.rbfe_utils import SolventComponentUtils
from openff.units import unit  # type: ignore
import pytest  # type: ignore
from openfe import SolventComponent  # type: ignore


def test_solvent_utils_loads_NoneType_input():
    with pytest.raises(TypeError):
        SolventComponentUtils.loads(None)


def test_solvent_utils_loads_invalid_input():
    # Invalid JSON string. The correct key is "positive_ion", and not "positive_io".
    invalid_keys_json = '{"positive_io": "Na", "negative_ion": "Cl", "neutralize": true, "ion_concentration": "0.15 mol/l"}'
    with pytest.raises(KeyError):
        SolventComponentUtils.loads(invalid_keys_json)

    # Invalid JSON string. NaN is not a valid value for positive_ion
    invalid_values_json = '{"positive_ion": "NaN", "negative_ion": "ClA", "neutralize": true, "ion_concentration": "0.15 mol/l"}'
    with pytest.raises(ValueError):
        SolventComponentUtils.loads(invalid_values_json)


def test_solvent_utils_loads_valid_input():

    # A valid solvent JSON string with all necessary keys.
    valid_json = '{"positive_ion": "Na", "negative_ion": "Cl", "neutralize": true, "ion_concentration": "0.15 mol/l"}'

    # Loading the JSON
    solvent = SolventComponentUtils.loads(valid_json)

    assert solvent.ion_concentration == 0.15 * unit.mole / unit.litre
    assert solvent.positive_ion == 'Na+'
    assert solvent.negative_ion == 'Cl-'
    assert solvent.neutralize is True


def test_solvent_utils_dumps_NoneType_input():
    with pytest.raises(TypeError):
        SolventComponentUtils.dumps(None)


def test_solvent_utils_dumps_invalid_input():
    # Testing with a dict object
    incorrect_type_object = {}

    with pytest.raises(TypeError):
        SolventComponentUtils.dumps(incorrect_type_object)


def test_solvent_utils_dumps_valid_input():

    # Creating a SolventComponent object
    solvent = SolventComponent(positive_ion='Na',
                               negative_ion='Cl',
                               neutralize=True,
                               ion_concentration=0.15 * unit.mole / unit.litre)

    # Serializing the SolventComponent object
    solvent_json = SolventComponentUtils.dumps(solvent)
    assert solvent_json == '{"positive_ion": "Na+", "negative_ion": "Cl-", "ion_concentration": "0.15 mole / liter", "neutralize": true}'
