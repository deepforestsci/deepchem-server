import os
from pathlib import Path


def get_cache_root() -> Path:
    """Return the cache root directory from the environment.

    Returns
    -------
    Path
        The cache root directory.

    Raises
    ------
    EnvironmentError
        If PROTEOMESCAN_CACHE_ROOT is not set.
    """
    root = os.environ.get("PROTEOMESCAN_CACHE_ROOT")
    if not root:
        raise EnvironmentError("PROTEOMESCAN_CACHE_ROOT environment variable is not set.")
    return Path(root)


def gene_root(scan_id: str, gene_name: str) -> Path:
    """Return the cache directory for a specific gene within a scan.

    Parameters
    ----------
    scan_id : str
        Scan run identifier.
    gene_name : str
        Gene name (e.g. 'MAP2K1').

    Returns
    -------
    Path
        Per-gene cache directory.
    """
    return get_cache_root() / scan_id / gene_name


def pdb_raw_path(scan_id: str, gene_name: str, pdb_id: str) -> Path:
    """Return the local path for a raw (unprocessed) PDB file.

    Parameters
    ----------
    scan_id : str
        Scan run identifier.
    gene_name : str
        Gene name.
    pdb_id : str
        PDB identifier (uppercase, eg '7M0X').

    Returns
    -------
    Path
        Expected local path of the raw PDB file.
    """
    return gene_root(scan_id, gene_name) / f"g_{gene_name}_p_{pdb_id}.pdb"


def pdb_clean_path(scan_id: str, gene_name: str, pdb_id: str) -> Path:
    """Return the local path for a cleaned PDB file.

    Parameters
    ----------
    scan_id : str
        Scan run identifier.
    gene_name : str
        Gene name.
    pdb_id : str
        PDB identifier (uppercase, eg '7M0X').

    Returns
    -------
    Path
        Expected local path of the cleaned PDB file.
    """
    return gene_root(scan_id, gene_name) / f"cleaned_g_{gene_name}_p_{pdb_id}.pdb"
