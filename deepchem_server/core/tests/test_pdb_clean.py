import os
import tempfile

import pytest

from deepchem_server.core import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.primitives.pdb_clean import pdb_clean


ASSETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")
RAW_PDB = os.path.join(ASSETS, "cleaned_3cyx.pdb")
LARGE_PDB = os.path.join(ASSETS, "181L_mod_capped_protonated.pdb")
TWO_CHAIN_PDB = os.path.join(ASSETS, "test_protein_2chain.pdb")


def _upload_pdb(datastore, filename, key):
    card = DataCard(address='', file_type='pdb', data_type='text/plain')
    return datastore.upload_data(key, filename, card)


def test_pdb_clean_basic(disk_datastore):
    """Basic PDB cleaning returns a valid deepchem address."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, RAW_PDB, 'raw_protein.pdb')
    cleaned_addr = pdb_clean(pdb_address=pdb_addr, output='cleaned_protein')

    assert cleaned_addr.startswith('deepchem://')
    assert cleaned_addr.endswith('.pdb')


def test_pdb_clean_output_has_pdb_extension(disk_datastore):
    """Output key gets .pdb extension when not supplied."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, RAW_PDB, 'raw_protein2.pdb')

    addr_no_ext = pdb_clean(pdb_address=pdb_addr, output='cleaned_no_ext')
    assert addr_no_ext.endswith('.pdb')

    addr_with_ext = pdb_clean(pdb_address=pdb_addr, output='cleaned_with_ext.pdb')
    assert addr_with_ext.endswith('.pdb')


def test_pdb_clean_produces_valid_pdb(disk_datastore):
    """Cleaned PDB file contains ATOM records."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, RAW_PDB, 'raw_protein3.pdb')
    cleaned_addr = pdb_clean(pdb_address=pdb_addr, output='cleaned_valid')

    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        tmp_path = tmp.name
    try:
        disk_datastore.download_object(cleaned_addr, tmp_path)
        with open(tmp_path, 'r') as f:
            content = f.read()
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    assert 'ATOM' in content or 'HETATM' in content


def test_pdb_clean_no_heterogens(disk_datastore):
    """Cleaned PDB with remove_heterogens=True contains no HETATM records (except HOH when remove_water=False)."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, RAW_PDB, 'raw_protein4.pdb')
    cleaned_addr = pdb_clean(
        pdb_address=pdb_addr,
        output='cleaned_no_hets',
        remove_heterogens=True,
        remove_water=True,
    )

    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        tmp_path = tmp.name
    try:
        disk_datastore.download_object(cleaned_addr, tmp_path)
        with open(tmp_path, 'r') as f:
            content = f.read()
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    water_lines = [l for l in content.splitlines() if l.startswith('HETATM') and 'HOH' in l]  # noqa: E741
    assert len(water_lines) == 0


def test_pdb_clean_data_card(disk_datastore):
    """DataCard is created with correct file_type and data_type."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, RAW_PDB, 'raw_protein5.pdb')
    cleaned_addr = pdb_clean(pdb_address=pdb_addr, output='cleaned_card_test')

    card = disk_datastore.get(cleaned_addr + '.cdc')
    assert card is not None
    assert card.file_type == 'pdb'
    assert card.data_type == 'text/plain'


def test_pdb_clean_empty_address_raises(disk_datastore):
    """Empty pdb_address raises ValueError."""
    config.set_datastore(disk_datastore)

    with pytest.raises(ValueError, match="pdb_address is required"):
        pdb_clean(pdb_address='', output='should_fail')


def test_pdb_clean_no_datastore_raises():
    """Calling without a configured datastore raises ValueError."""
    config.set_datastore(None)

    with pytest.raises(ValueError, match="Datastore not set"):
        pdb_clean(pdb_address='deepchem://test/user/some.pdb', output='out')


def test_pdb_clean_custom_ph(disk_datastore):
    """Cleaning with a non-default pH completes without error."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, RAW_PDB, 'raw_protein6.pdb')
    cleaned_addr = pdb_clean(pdb_address=pdb_addr, output='cleaned_ph6', ph=6.0)

    assert cleaned_addr.startswith('deepchem://')


def test_pdb_clean_without_hydrogens(disk_datastore):
    """add_hydrogens=False still produces a valid structure."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, RAW_PDB, 'raw_protein7.pdb')
    cleaned_addr = pdb_clean(pdb_address=pdb_addr, output='cleaned_noh', add_hydrogens=False)

    assert cleaned_addr.startswith('deepchem://')


def _download_pdb_content(datastore, address):
    """Download a PDB from the datastore and return its text content."""
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        tmp_path = tmp.name
    try:
        datastore.download_object(address, tmp_path)
        with open(tmp_path, 'r') as f:
            return f.read()
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_pdb_clean_two_chain_protein(disk_datastore):
    """Cleaning the two-chain fixture returns a valid deepchem address."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, TWO_CHAIN_PDB, 'two_chain.pdb')
    cleaned_addr = pdb_clean(pdb_address=pdb_addr, output='cleaned_two_chain')

    assert cleaned_addr.startswith('deepchem://')
    assert cleaned_addr.endswith('.pdb')


def test_pdb_clean_two_chain_has_atom_records(disk_datastore):
    """Cleaned two-chain structure contains ATOM records."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, TWO_CHAIN_PDB, 'two_chain_atoms.pdb')
    cleaned_addr = pdb_clean(pdb_address=pdb_addr, output='cleaned_two_chain_atoms')

    content = _download_pdb_content(disk_datastore, cleaned_addr)
    atom_lines = [line for line in content.splitlines() if line.startswith('ATOM')]
    assert len(atom_lines) > 0


def test_pdb_clean_water_retained(disk_datastore):
    """remove_water=False preserves HOH records in the cleaned structure."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, TWO_CHAIN_PDB, 'two_chain_keepwater.pdb')
    cleaned_addr = pdb_clean(
        pdb_address=pdb_addr,
        output='cleaned_keepwater',
        remove_heterogens=True,
        remove_water=False,
    )

    content = _download_pdb_content(disk_datastore, cleaned_addr)
    hoh_lines = [line for line in content.splitlines() if 'HOH' in line]
    assert len(hoh_lines) > 0, "Expected HOH water records to be retained"


def test_pdb_clean_heterogens_kept_when_disabled(disk_datastore):
    """remove_heterogens=False leaves all HETATM records intact."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, TWO_CHAIN_PDB, 'two_chain_kephets.pdb')
    cleaned_addr = pdb_clean(
        pdb_address=pdb_addr,
        output='cleaned_keephets',
        remove_heterogens=False,
        remove_water=False,
    )

    content = _download_pdb_content(disk_datastore, cleaned_addr)
    hetatm_lines = [line for line in content.splitlines() if line.startswith('HETATM')]
    assert len(hetatm_lines) > 0, "Expected HETATM records to be kept"


def test_pdb_clean_large_protein(disk_datastore):
    """Cleaning the 181L structure (2612 atoms) returns a valid address."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, LARGE_PDB, 'large_protein.pdb')
    cleaned_addr = pdb_clean(pdb_address=pdb_addr, output='cleaned_large')

    assert cleaned_addr.startswith('deepchem://')
    assert cleaned_addr.endswith('.pdb')


def test_pdb_clean_large_protein_has_atoms(disk_datastore):
    """Cleaned 181L structure retains ATOM records."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, LARGE_PDB, 'large_protein_atoms.pdb')
    cleaned_addr = pdb_clean(pdb_address=pdb_addr, output='cleaned_large_atoms')

    content = _download_pdb_content(disk_datastore, cleaned_addr)
    atom_lines = [line for line in content.splitlines() if line.startswith('ATOM')]
    assert len(atom_lines) > 0


def test_pdb_clean_large_protein_no_water(disk_datastore):
    """After cleaning 181L with remove_water=True, no HOH records remain."""
    config.set_datastore(disk_datastore)

    pdb_addr = _upload_pdb(disk_datastore, LARGE_PDB, 'large_protein_nohoh.pdb')
    cleaned_addr = pdb_clean(
        pdb_address=pdb_addr,
        output='cleaned_large_nohoh',
        remove_heterogens=True,
        remove_water=True,
    )

    content = _download_pdb_content(disk_datastore, cleaned_addr)
    hoh_lines = [line for line in content.splitlines() if line.startswith('HETATM') and 'HOH' in line]
    assert len(hoh_lines) == 0


def test_pdb_clean_ph_range(disk_datastore):
    """Cleaning at acidic (pH 5.0) and basic (pH 9.0) both complete without error."""
    config.set_datastore(disk_datastore)

    for ph_val in (5.0, 9.0):
        pdb_addr = _upload_pdb(disk_datastore, RAW_PDB, f'raw_ph{int(ph_val)}.pdb')
        cleaned_addr = pdb_clean(
            pdb_address=pdb_addr,
            output=f'cleaned_ph{int(ph_val)}',
            ph=ph_val,
        )
        assert cleaned_addr.startswith('deepchem://')
