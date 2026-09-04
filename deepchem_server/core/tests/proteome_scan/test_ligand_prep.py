import json
import os

import pytest
from rdkit import Chem

from deepchem_server.core import config
from deepchem_server.core.primitives.proteome_scan.ligand_prep import ligand_prep


ASSETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets")

with open(os.path.join(ASSETS, "ligand_SMILES.json")) as _f:
    _LIGAND_DATA = json.load(_f)

SMILES_1 = 'CC(=O)Oc1ccccc1C(=O)O'
SMILES_2 = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O'
SMILES_3 = 'CCO'

PROTEOME_SCAN_SUBSET = {
    "Erlotinib": _LIGAND_DATA["ligand_smiles"]["Erlotinib"],
    "Selinexor": _LIGAND_DATA["ligand_smiles"]["Selinexor"],
    "Palbociclib": _LIGAND_DATA["ligand_smiles"]["Palbociclib"],
    "Olaparib": _LIGAND_DATA["ligand_smiles"]["Olaparib"],
}

CHIRAL_NAME = "SOS1-IN-11"
CHIRAL_SMILES = _LIGAND_DATA["ligand_smiles"][CHIRAL_NAME]


def test_ligand_prep_basic(disk_datastore):
    """Basic SDF generation from a valid SMILES."""
    config.set_datastore(disk_datastore)

    address = ligand_prep(smiles=SMILES_1, output='smiles_1')

    assert address.startswith('deepchem://')
    assert address.endswith('.sdf')


def test_ligand_prep_output_has_sdf_extension(disk_datastore):
    """Output key gets .sdf extension when not provided."""
    config.set_datastore(disk_datastore)

    address = ligand_prep(smiles=SMILES_3, output='smiles_3')
    assert address.endswith('.sdf')

    address2 = ligand_prep(smiles=SMILES_3, output='smiles_32.sdf')
    assert address2.endswith('.sdf')


def test_ligand_prep_sdf_content_is_valid(disk_datastore):
    """Retrieved SDF file contains expected molecule keywords."""
    config.set_datastore(disk_datastore)

    address = ligand_prep(smiles=SMILES_1, output='smiles_1_content_test')

    sdf_lines = disk_datastore.get(address)
    sdf_text = ''.join(sdf_lines) if isinstance(sdf_lines, list) else sdf_lines
    assert 'M  END' in sdf_text
    assert '$$$$' in sdf_text


def test_ligand_prep_with_ligand_name(disk_datastore):
    """Ligand name is embedded in the SDF block."""
    config.set_datastore(disk_datastore)

    address = ligand_prep(smiles=SMILES_1, output='smiles_1_named', ligand_name='smiles_1')

    sdf_lines = disk_datastore.get(address)
    sdf_text = ''.join(sdf_lines) if isinstance(sdf_lines, list) else sdf_lines
    assert 'smiles_1' in sdf_text


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


def test_ligand_prep_proteome_scan_drugs(disk_datastore):
    """Subset of 4 ProteomeScan clinical drug molecules all produce valid SDF."""
    config.set_datastore(disk_datastore)

    for name, smiles in PROTEOME_SCAN_SUBSET.items():
        address = ligand_prep(smiles=smiles, output=f'ps_{name.lower()}', ligand_name=name)
        assert address.startswith('deepchem://')
        assert address.endswith('.sdf')

        sdf_lines = disk_datastore.get(address)
        sdf_text = ''.join(sdf_lines) if isinstance(sdf_lines, list) else sdf_lines
        assert 'M  END' in sdf_text, f"M  END missing from SDF for {name}"
        assert '$$$$' in sdf_text, f"$$$$ record missing from SDF for {name}"


def test_ligand_prep_all_dataset_smiles_parseable():
    """Every non-empty SMILES in the ProteomeScan dataset is valid RDKit-parseable."""
    failures = []
    for name, smiles in _LIGAND_DATA["ligand_smiles"].items():
        if not smiles:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            failures.append(name)
    assert failures == [], f"Invalid SMILES in dataset for: {failures}"
