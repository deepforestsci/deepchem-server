"""
Server-side cache root for proteome scan artifacts.
Directory layout under the cache root:

    <cache_root>/<scan_id>/<gene_name>/g_<gene>_p_<pdb>.pdb
    <cache_root>/<scan_id>/<gene_name>/cleaned_g_<gene>_p_<pdb>.pdb
    <cache_root>/<scan_id>/<gene_name>/<gene>_pdbs.csv
    <cache_root>/<scan_id>/<gene_name>/complexes/complex_<gene>_<pdb>_<ligand>.pdb
    <cache_root>/<scan_id>/<ligand>/top_score_<gene>_<ligand>.csv
    <cache_root>/<scan_id>/scan_results/top_score_<ligand>.csv

The cache root is resolved in this order:

1. PROTEOMESCAN_CACHE_ROOT
2. $DEEPCHEM_CACHE_HOME/deepchem/proteome_scan
3. ~/.cache/deepchem/proteome_scan
"""

import os
import shutil
from pathlib import Path
from typing import Optional


DEFAULT_SUBDIR = "proteome_scan"


def get_cache_root() -> Path:
    """
    Get the base cache directory for all proteome scan artifacts.

    The location is resolved in this order:
    1) PROTEOMESCAN_CACHE_ROOT
    2) DEEPCHEM_CACHE_HOME/deepchem/proteome_scan
    3) ~/.cache/deepchem/proteome_scan

    Returns
    -------
    pathlib.Path
        Cache root directory. The directory is created if needed.

    Examples
    --------
    >>> root = get_cache_root()
    >>> root.exists()
    True
    """
    override = os.getenv("PROTEOMESCAN_CACHE_ROOT")
    if override:
        root = Path(override).expanduser()
    else:
        cache_home = os.getenv("DEEPCHEM_CACHE_HOME")
        if cache_home:
            root = Path(cache_home).expanduser() / "deepchem" / DEFAULT_SUBDIR
        else:
            root = Path.home() / ".cache" / "deepchem" / DEFAULT_SUBDIR
    root.mkdir(parents=True, exist_ok=True)
    return root


def get_scan_dir(scan_id: str) -> Path:
    """
    Get the on-disk directory for a scan.

    Parameters
    ----------
    scan_id : str
        Scan identifier.

    Returns
    -------
    pathlib.Path
        Scan directory. The directory is created if needed.

    Examples
    --------
    >>> d = get_scan_dir("scan1")
    >>> d.name
    'scan1'
    """
    scan_dir = get_cache_root() / scan_id
    scan_dir.mkdir(parents=True, exist_ok=True)
    return scan_dir


def get_gene_dir(scan_id: str, gene_name: str) -> Path:
    """
    Get the per-gene directory inside a scan.

    Parameters
    ----------
    scan_id : str
        Scan identifier.
    gene_name : str
        Gene symbol.

    Returns
    -------
    pathlib.Path
        Per-gene directory. The directory is created if needed.

    Examples
    --------
    >>> d = get_gene_dir("scan1", "TP53")
    >>> d.name
    'TP53'
    """
    gene_dir = get_scan_dir(scan_id) / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)
    return gene_dir


def get_gene_complexes_dir(scan_id: str, gene_name: str) -> Path:
    """
    Get the directory where per-gene complex PDBs are stored.

    Parameters
    ----------
    scan_id : str
        Scan identifier.
    gene_name : str
        Gene symbol.

    Returns
    -------
    pathlib.Path
        Complexes directory. The directory is created if needed.

    Examples
    --------
    >>> d = get_gene_complexes_dir("scan1", "TP53")
    >>> d.name
    'complexes'
    """
    complexes_dir = get_gene_dir(scan_id, gene_name) / "complexes"
    complexes_dir.mkdir(parents=True, exist_ok=True)
    return complexes_dir


def get_ligand_dir(scan_id: str, ligand_name: str) -> Path:
    """
    Get the per-ligand directory inside a scan.

    Parameters
    ----------
    scan_id : str
        Scan identifier.
    ligand_name : str
        Ligand label.

    Returns
    -------
    pathlib.Path
        Per-ligand directory. The directory is created if needed.

    Examples
    --------
    >>> d = get_ligand_dir("scan1", "LIG")
    >>> d.name
    'LIG'
    """
    ligand_dir = get_scan_dir(scan_id) / ligand_name
    ligand_dir.mkdir(parents=True, exist_ok=True)
    return ligand_dir


def get_scan_results_dir(scan_id: str) -> Path:
    """
    Get the directory for aggregated scan-level outputs.

    Parameters
    ----------
    scan_id : str
        Scan identifier.

    Returns
    -------
    pathlib.Path
        Directory for scan-level results. The directory is created if needed.

    Examples
    --------
    >>> d = get_scan_results_dir("scan1")
    >>> d.name
    'scan_results'
    """
    results_dir = get_scan_dir(scan_id) / "scan_results"
    results_dir.mkdir(parents=True, exist_ok=True)
    return results_dir


def raw_pdb_path(scan_id: str, gene_name: str, pdb_id: str) -> Path:
    """
    Get the expected path for a raw downloaded PDB file.

    Parameters
    ----------
    scan_id : str
        Scan identifier.
    gene_name : str
        Gene symbol.
    pdb_id : str
        PDB identifier.

    Returns
    -------
    pathlib.Path
        Path under the per-gene scan directory.

    Examples
    --------
    >>> p = raw_pdb_path("scan1", "TP53", "1ABC")
    >>> p.name.endswith(".pdb")
    True
    """
    return get_gene_dir(scan_id, gene_name) / f"g_{gene_name}_p_{pdb_id}.pdb"


def cleaned_pdb_path(scan_id: str, gene_name: str, pdb_id: str) -> Path:
    """
    Get the expected path for a cleaned PDB file.

    Parameters
    ----------
    scan_id : str
        Scan identifier.
    gene_name : str
        Gene symbol.
    pdb_id : str
        PDB identifier.

    Returns
    -------
    pathlib.Path
        Path under the per-gene scan directory.

    Examples
    --------
    >>> p = cleaned_pdb_path("scan1", "TP53", "1ABC")
    >>> "cleaned" in p.name
    True
    """
    return get_gene_dir(scan_id, gene_name) / f"cleaned_g_{gene_name}_p_{pdb_id}.pdb"


def gene_pdbs_csv_path(scan_id: str, gene_name: str) -> Path:
    """
    Get the expected path for the per-gene PDB metadata CSV.

    Parameters
    ----------
    scan_id : str
        Scan identifier.
    gene_name : str
        Gene symbol.

    Returns
    -------
    pathlib.Path
        CSV path under the per-gene scan directory.

    Examples
    --------
    >>> p = gene_pdbs_csv_path("scan1", "TP53")
    >>> p.suffix
    '.csv'
    """
    return get_gene_dir(scan_id, gene_name) / f"{gene_name}_pdbs.csv"


def complex_path(scan_id: str, gene_name: str, pdb_id: str, ligand_name: str) -> Path:
    """
    Get the expected path for a docked complex PDB file.

    Parameters
    ----------
    scan_id : str
        Scan identifier.
    gene_name : str
        Gene symbol.
    pdb_id : str
        PDB identifier.
    ligand_name : str
        Ligand label.

    Returns
    -------
    pathlib.Path
        Path under the per-gene complexes directory.

    Examples
    --------
    >>> p = complex_path("scan1", "TP53", "1ABC", "LIG")
    >>> p.name.startswith("complex_")
    True
    """
    return (get_gene_complexes_dir(scan_id, gene_name) / f"complex_{gene_name}_{pdb_id}_{ligand_name}.pdb")


def top_score_gene_ligand_csv_path(scan_id: str, gene_name: str, ligand_name: str) -> Path:
    """
    Get the expected path for the per-gene, per-ligand score CSV.

    Parameters
    ----------
    scan_id : str
        Scan identifier.
    gene_name : str
        Gene symbol.
    ligand_name : str
        Ligand label.

    Returns
    -------
    pathlib.Path
        CSV path under the per-ligand scan directory.

    Examples
    --------
    >>> p = top_score_gene_ligand_csv_path("scan1", "TP53", "LIG")
    >>> p.suffix
    '.csv'
    """
    return get_ligand_dir(scan_id, ligand_name) / f"top_score_{gene_name}_{ligand_name}.csv"


def top_score_ligand_csv_path(scan_id: str, ligand_name: str) -> Path:
    """
    Get the expected path for the aggregated per-ligand score CSV.

    Parameters
    ----------
    scan_id : str
        Scan identifier.
    ligand_name : str
        Ligand label.

    Returns
    -------
    pathlib.Path
        CSV path under the scan results directory.

    Examples
    --------
    >>> p = top_score_ligand_csv_path("scan1", "LIG")
    >>> p.name.startswith("top_score_")
    True
    """
    return get_scan_results_dir(scan_id) / f"top_score_{ligand_name}.csv"


def clear_scan(scan_id: str) -> None:
    """
    Delete a scan directory and everything under it.

    Parameters
    ----------
    scan_id : str
        Scan identifier to delete.

    Returns
    -------
    None

    Examples
    --------
    >>> _ = get_scan_dir("scan_to_delete")
    >>> clear_scan("scan_to_delete")
    """
    scan_dir = get_cache_root() / scan_id
    if scan_dir.exists():
        shutil.rmtree(scan_dir)


def list_gene_ligand_csvs(scan_id: str, ligand_name: str) -> list:
    """
    List per-gene score CSV files for a ligand.

    Parameters
    ----------
    scan_id : str
        Scan identifier.
    ligand_name : str
        Ligand label.

    Returns
    -------
    list
        Sorted list of file paths.

    Examples
    --------
    >>> paths = list_gene_ligand_csvs("scan1", "LIG")
    >>> isinstance(paths, list)
    True
    """
    ligand_dir = get_ligand_dir(scan_id, ligand_name)
    return sorted(str(p) for p in ligand_dir.iterdir() if p.is_file() and p.suffix == ".csv")


def ensure_parent(path: Path) -> Path:
    """
    Ensure the parent directory of a path exists.

    Parameters
    ----------
    path : pathlib.Path
        Path whose parent directory should exist.

    Returns
    -------
    pathlib.Path
        The input path.

    Examples
    --------
    >>> p = ensure_parent(Path("/tmp/deepchem_server/a/b.txt"))
    >>> p.parent.exists()
    True
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def safe_mkdir(path: Path, exist_ok: bool = True) -> Path:
    """
    Create a directory (and parents) if needed.

    Parameters
    ----------
    path : pathlib.Path
        Directory path to create.
    exist_ok : bool, default True
        If True, do not error when the directory already exists.

    Returns
    -------
    pathlib.Path
        The created directory path.

    Examples
    --------
    >>> d = safe_mkdir(Path("/tmp/deepchem_server/cache_test"))
    >>> d.exists()
    True
    """
    path.mkdir(parents=True, exist_ok=exist_ok)
    return path


def get_optional_env(name: str) -> Optional[str]:
    """Helper used by tests to introspect the active cache env var."""
    return os.getenv(name)
