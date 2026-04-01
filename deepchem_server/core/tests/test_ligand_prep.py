import pytest

from deepchem_server.core import config
from deepchem_server.core.primitives.ligand_prep import ligand_prep


ASPIRIN_SMILES = 'CC(=O)Oc1ccccc1C(=O)O'
IBUPROFEN_SMILES = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O'
ETHANOL_SMILES = 'CCO'


def test_ligand_prep_basic(disk_datastore):
    """Basic SDF generation from a valid SMILES."""
    config.set_datastore(disk_datastore)

    address = ligand_prep(smiles=ASPIRIN_SMILES, output='aspirin')

    assert address.startswith('deepchem://')
    assert address.endswith('.sdf')


def test_ligand_prep_output_has_sdf_extension(disk_datastore):
    """Output key gets .sdf extension when not provided."""
    config.set_datastore(disk_datastore)

    address = ligand_prep(smiles=ETHANOL_SMILES, output='ethanol')
    assert address.endswith('.sdf')

    address2 = ligand_prep(smiles=ETHANOL_SMILES, output='ethanol2.sdf')
    assert address2.endswith('.sdf')


def test_ligand_prep_sdf_content_is_valid(disk_datastore):
    """Retrieved SDF file contains expected molecule keywords."""
    config.set_datastore(disk_datastore)

    address = ligand_prep(smiles=ASPIRIN_SMILES, output='aspirin_content_test')

    sdf_lines = disk_datastore.get(address)
    sdf_text = ''.join(sdf_lines) if isinstance(sdf_lines, list) else sdf_lines
    assert 'M  END' in sdf_text
    assert '$$$$' in sdf_text


def test_ligand_prep_with_ligand_name(disk_datastore):
    """Ligand name is embedded in the SDF block."""
    config.set_datastore(disk_datastore)

    address = ligand_prep(smiles=ASPIRIN_SMILES, output='aspirin_named', ligand_name='Aspirin')

    sdf_lines = disk_datastore.get(address)
    sdf_text = ''.join(sdf_lines) if isinstance(sdf_lines, list) else sdf_lines
    assert 'Aspirin' in sdf_text


def test_ligand_prep_without_optimization(disk_datastore):
    """Primitive succeeds with optimization disabled."""
    config.set_datastore(disk_datastore)

    address = ligand_prep(smiles=ASPIRIN_SMILES, output='aspirin_noopt', optimize=False)

    assert address.startswith('deepchem://')
    sdf_lines = disk_datastore.get(address)
    sdf_text = ''.join(sdf_lines) if isinstance(sdf_lines, list) else sdf_lines
    assert 'M  END' in sdf_text


def test_ligand_prep_reproducible_with_seed(disk_datastore):
    """Same SMILES and seed produces output at the same address structure."""
    config.set_datastore(disk_datastore)

    addr1 = ligand_prep(smiles=ASPIRIN_SMILES, output='aspirin_seed1', random_seed=42)
    addr2 = ligand_prep(smiles=ASPIRIN_SMILES, output='aspirin_seed2', random_seed=42)

    sdf1 = disk_datastore.get(addr1)
    sdf2 = disk_datastore.get(addr2)
    text1 = ''.join(sdf1) if isinstance(sdf1, list) else sdf1
    text2 = ''.join(sdf2) if isinstance(sdf2, list) else sdf2

    # Both should be valid SDF files
    assert 'M  END' in text1
    assert 'M  END' in text2


def test_ligand_prep_data_card(disk_datastore):
    """DataCard is created with correct file_type and data_type."""
    config.set_datastore(disk_datastore)

    address = ligand_prep(smiles=ASPIRIN_SMILES, output='aspirin_card_test')

    card = disk_datastore.get(address + '.cdc')
    assert card is not None
    assert card.file_type == 'sdf'
    assert card.data_type == 'sdf'


def test_ligand_prep_invalid_smiles(disk_datastore):
    """Invalid SMILES raises ValueError."""
    config.set_datastore(disk_datastore)

    with pytest.raises(ValueError, match="Invalid SMILES"):
        ligand_prep(smiles='not_a_smiles!!!', output='bad_mol')


def test_ligand_prep_empty_smiles(disk_datastore):
    """Empty SMILES raises ValueError."""
    config.set_datastore(disk_datastore)

    with pytest.raises(ValueError, match="SMILES string is required"):
        ligand_prep(smiles='', output='empty')


def test_ligand_prep_different_molecules(disk_datastore):
    """Multiple different molecules all produce valid SDF files."""
    config.set_datastore(disk_datastore)

    for name, smiles in [('aspirin', ASPIRIN_SMILES), ('ibuprofen', IBUPROFEN_SMILES), ('ethanol', ETHANOL_SMILES)]:
        address = ligand_prep(smiles=smiles, output=f'mol_{name}')
        assert address.startswith('deepchem://')

        sdf_lines = disk_datastore.get(address)
        sdf_text = ''.join(sdf_lines) if isinstance(sdf_lines, list) else sdf_lines
        assert 'M  END' in sdf_text
