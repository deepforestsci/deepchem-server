import os
import tempfile
from typing import List, Optional

from deepchem_server.core.common import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.common.progress_logger import log_progress


def pdb_clean(
    pdb_address: str,
    output: str,
    remove_chains: Optional[List[str]] = None,
    remove_heterogens: bool = True,
    remove_water: bool = True,
    add_hydrogens: bool = True,
    ph: float = 7.0,
) -> str:
    """
    Clean a PDB file using PDBFixer and OpenMM.

    Removes heterogens and water, optionally removes extra chains, and
    adds missing hydrogens at a given pH.

    1. Load the raw PDB via PDBFixer.
    2. (Optional) Remove specified chains.
    3. Remove heterogens if remove_heterogens is True (with or without water) (default True).
    4. Add missing hydrogens if add_hydrogens is True at ph (default True).
    5. Write the cleaned structure with OpenMM PDBFile.

    Parameters
    ----------
    pdb_address : str
        DeepChem address of the input PDB file.
    output : str
        Output filename key in the datastore (with or without .pdb extension).
    remove_chains : list of str, optional
        Chain IDs to remove from the structure before cleaning.
    remove_heterogens : bool, optional
        Whether to strip heteroatom records (default True).
    remove_water : bool, optional
        Whether to also remove water molecules (default True).
        Only has effect when remove_heterogens is True.
    add_hydrogens : bool, optional
        Whether to add missing hydrogens (default True).
    ph : float, optional
        pH for hydrogen protonation (default 7.0).

    Returns
    -------
    str
        DeepChem address of the cleaned PDB file.

    Raises
    ------
    ImportError
        If PDBFixer or OpenMM are not installed.
    ValueError
        If pdb_address is empty or the datastore is not set.
    RuntimeError
        If PDBFixer fails to clean the structure.

    Examples
    --------
    >>> cleaned_addr = pdb_clean(
    ...     pdb_address="deepchem://profile/project/raw_protein.pdb",
    ...     output="cleaned_protein",
    ... )
    >>> print(cleaned_addr)
    deepchem://profile/project/cleaned_protein.pdb
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")

    if not pdb_address:
        raise ValueError("pdb_address is required")

    try:
        from pdbfixer.pdbfixer import PDBFixer  # type: ignore
    except ImportError:
        raise ImportError("pdbfixer is required for pdb_clean but not installed")

    try:
        from openmm.app import PDBFile  # type: ignore
    except ImportError:
        raise ImportError("OpenMM is required for pdb_clean but not installed")

    output_key = output if output.endswith('.pdb') else f'{output}.pdb'

    with tempfile.TemporaryDirectory() as tmp:
        log_progress('pdb_clean', 10, f'downloading PDB from {pdb_address}')
        raw_pdb_path = os.path.join(tmp, 'input.pdb')
        datastore.download_object(pdb_address, raw_pdb_path)

        log_progress('pdb_clean', 30, 'loading structure with PDBFixer')
        fixer = PDBFixer(raw_pdb_path)

        if remove_chains:
            log_progress('pdb_clean', 40, f'removing chains: {remove_chains}')
            fixer.removeChains(chainIds=remove_chains)

        if remove_heterogens:
            keep_water = not remove_water
            log_progress('pdb_clean', 50, f'removing heterogens (keep_water={keep_water})')
            fixer.removeHeterogens(keep_water)

        if add_hydrogens:
            log_progress('pdb_clean', 65, f'adding missing hydrogens at pH {ph}')
            try:
                fixer.addMissingHydrogens(ph)
            except Exception:
                fixer.findMissingResidues()
                fixer.findMissingAtoms()
                fixer.addMissingAtoms()
                fixer.addMissingHydrogens(ph)

        log_progress('pdb_clean', 80, 'writing cleaned PDB')
        cleaned_pdb_path = os.path.join(tmp, 'cleaned.pdb')
        try:
            with open(cleaned_pdb_path, 'w') as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)
        except Exception as e:
            raise RuntimeError(f"Failed to write cleaned structure: {e}")

        with open(cleaned_pdb_path, 'r') as f:
            cleaned_content = f.read()

    log_progress('pdb_clean', 92, 'uploading cleaned PDB to datastore')
    card = DataCard(address='', file_type='pdb', data_type='text/plain')
    address = datastore.upload_data_from_memory(cleaned_content, output_key, card)

    log_progress('pdb_clean', 100, 'PDB cleaning complete')
    return address
