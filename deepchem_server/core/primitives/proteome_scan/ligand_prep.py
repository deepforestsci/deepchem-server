"""
Prepare a ligand for docking.
"""

import io
from typing import Optional

from deepchem_server.core.common import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.common.progress_logger import log_progress

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


def ligand_prep(
    smiles: str,
    output: str,
    ligand_name: Optional[str] = "",
    random_seed: Optional[int] = None,
) -> str:
    """
    Convert a SMILES string to a 3D SDF file using RDKit.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    output : str
        Output filename in the datastore (with or without .sdf extension).
    ligand_name : str, optional
        Name to embed in the SDF molecule block.
    random_seed : int, optional
        Reproducibility seed for 3D embedding.

    Returns
    -------
    str
        DeepChem address of the generated SDF file.

    Examples
    --------
    >>> addr = ligand_prep("CCO", output="ethanol", ligand_name="EtOH")
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")

    if not smiles:
        raise ValueError("SMILES string is required")

    if not RDKIT_AVAILABLE:
        raise ImportError(
            "RDKit is required for ligand preparation but not installed")

    log_progress('ligand_prep', 10, f'parsing SMILES: {smiles}')
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    if ligand_name:
        mol.SetProp('_Name', ligand_name)

    log_progress('ligand_prep', 30, 'adding hydrogens')
    mol = Chem.AddHs(mol)

    log_progress("ligand_prep", 50,
                 "generating 3D coordinates with ETKDG")
    params = AllChem.ETKDGv3()  # type: ignore
    if random_seed is not None:
        params.randomSeed = random_seed
    result = AllChem.EmbedMolecule(mol, params)  # type: ignore
    if result == -1:
        raise ValueError(f"3D embedding failed for SMILES: {smiles}")

    log_progress("ligand_prep", 80, "serializing to SDF")
    sdf_buffer = io.StringIO()
    writer = Chem.SDWriter(sdf_buffer)
    writer.write(mol)
    writer.close()
    sdf_content = sdf_buffer.getvalue()

    output_key = output if output.endswith('.sdf') else f'{output}.sdf'
    card = DataCard(address='', file_type='sdf', data_type='sdf')

    log_progress('ligand_prep', 95, 'uploading SDF to datastore')
    address = datastore.upload_data_from_memory(sdf_content, output_key,
                                                card)

    log_progress('ligand_prep', 100, 'ligand preparation complete')
    return address
