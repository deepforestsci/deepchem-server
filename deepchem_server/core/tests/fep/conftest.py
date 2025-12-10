from pathlib import Path

import pytest  # type: ignore

from deepchem_server.core import config
from deepchem_server.core.common.cards import DataCard


ASSETS_DIR = Path(__file__).resolve().parent.parent / "assets"


@pytest.fixture(autouse=True)
def set_datastore(disk_datastore):
    config.set_datastore(disk_datastore)


@pytest.fixture()
def ligands_datastore_address():
    LIGANDS_SDF_PATH = ASSETS_DIR / "ligands.sdf"
    ligands_data_card = DataCard(address='', file_type='sdf', data_type='text/plain')
    ligands_datastore_address = config.get_datastore().upload_data('ligands.sdf', LIGANDS_SDF_PATH.as_posix(),
                                                                   ligands_data_card)
    return ligands_datastore_address


@pytest.fixture()
def ligands_nested_datastore_address():
    LIGANDS_SDF_PATH = ASSETS_DIR / "ligands.sdf"
    ligands_data_card = DataCard(address='', file_type='sdf', data_type='text/plain')
    ligands_datastore_address = config.get_datastore().upload_data('test rbfe/ligands.sdf', LIGANDS_SDF_PATH.as_posix(),
                                                                   ligands_data_card)
    return ligands_datastore_address


@pytest.fixture()
def ligands_large_datastore_address():
    LIGANDS_SDF_PATH = ASSETS_DIR / "ligands_large.sdf"
    ligands_data_card = DataCard(address='', file_type='sdf', data_type='text/plain')
    ligands_datastore_address = config.get_datastore().upload_data('ligands.sdf', LIGANDS_SDF_PATH.as_posix(),
                                                                   ligands_data_card)
    return ligands_datastore_address


@pytest.fixture()
def ligands_nested_large_datastore_address():
    LIGANDS_SDF_PATH = ASSETS_DIR / "ligands_large.sdf"
    ligands_data_card = DataCard(address='', file_type='sdf', data_type='text/plain')
    ligands_datastore_address = config.get_datastore().upload_data('test rbfe/ligands.sdf', LIGANDS_SDF_PATH.as_posix(),
                                                                   ligands_data_card)
    return ligands_datastore_address


@pytest.fixture()
def protein_datastore_address():
    PROTEIN_SDF_PATH = ASSETS_DIR / "181L_mod_capped_protonated.pdb"
    protein_data_card = DataCard(address='', file_type='pdb', data_type='text/plain')
    protein_datastore_address = config.get_datastore().upload_data('cleaned-protein.pdb', PROTEIN_SDF_PATH.as_posix(),
                                                                   protein_data_card)
    return protein_datastore_address


@pytest.fixture()
def protein_nested_datastore_address():
    PROTEIN_SDF_PATH = ASSETS_DIR / "181L_mod_capped_protonated.pdb"
    protein_data_card = DataCard(address='', file_type='pdb', data_type='text/plain')
    protein_datastore_address = config.get_datastore().upload_data('test rbfe/cleaned-protein.pdb',
                                                                   PROTEIN_SDF_PATH.as_posix(), protein_data_card)
    return protein_datastore_address
