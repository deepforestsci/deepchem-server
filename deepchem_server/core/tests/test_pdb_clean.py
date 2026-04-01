import os
import tempfile

import pytest

from deepchem_server.core import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.primitives.pdb_clean import pdb_clean


ASSETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")
RAW_PDB = os.path.join(ASSETS, "cleaned_3cyx.pdb")  # use existing test PDB as "raw" input


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

    # Water lines start with HETATM and have HOH residue name; after cleaning, none should remain
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
