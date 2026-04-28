"""
Run docking for the proteome scan workflow.
"""

import datetime as dt
import json
import os
import tempfile
from typing import Any, Dict, List, Optional, Tuple, Union

from deepchem.dock.pose_generation import DOCKED_POSES
import pandas as pd

from deepchem_server.core.common import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.common.progress_logger import log_progress
from deepchem_server.core.primitives.proteome_scan import cache as ps_cache


def _vina_docking(
    raw_pdb_path: str,
    raw_ligand_path: str,
    work_dir: str,
    exhaustiveness: int = 32,
    num_modes: int = 8,
    **kwargs: Any,
) -> Union[Tuple[DOCKED_POSES, List[float]], DOCKED_POSES]:
    """
    Dock a ligand into a protein structure with Vina.

    This is a thin wrapper around DeepChem's VinaPoseGenerator.

    Parameters
    ----------
    raw_pdb_path : str
        Path to the protein PDB file.
    raw_ligand_path : str
        Path to the ligand SDF file.
    work_dir : str
        Working directory where intermediate outputs may be written.
    exhaustiveness : int, default 32
        Search exhaustiveness passed to the docking engine.
    num_modes : int, default 8
        Maximum number of poses to generate.
    **kwargs : Any
        Additional keyword arguments passed to be passed to deepchem.dock.generate_poses method.

    Returns
    -------
    tuple
        Tuple of (complex, scores) as returned by the pose generator.

    Examples
    --------
    >>> complex, scores = _vina_docking("protein.pdb", "ligand.sdf", "/tmp")
    """
    try:
        from deepchem.dock.pose_generation import (
            VinaPoseGenerator,)  # type: ignore
    except ImportError:
        raise ImportError("DeepChem is required for run_docking but not installed")

    dir_name = "vina_docking_" + str(dt.datetime.now().isoformat())
    tmp = os.path.join(work_dir, dir_name)
    os.makedirs(tmp, exist_ok=True)

    pg = VinaPoseGenerator()
    complex, scores = pg.generate_poses(
        molecular_complex=(raw_pdb_path, raw_ligand_path),
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        out_dir=tmp,
        generate_scores=True,
        **kwargs,
    )
    return complex, scores


def _download_ligand_to_cache(
    ligand_address: str,
    scan_id: str,
    ligand_name: str,
) -> str:
    """
    Download a ligand from the datastore into the local scan cache.

    Parameters
    ----------
    ligand_address : str
        Datastore address for the ligand SDF.
    scan_id : str
        Identifier for the scan cache namespace.
    ligand_name : str
        Ligand label used as a folder name and filename stem.

    Returns
    -------
    str
        Local path to the downloaded SDF file.

    Examples
    --------
    >>> path = _download_ligand_to_cache("ds://ligand.sdf", "scan1", "LIG")
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")

    ligand_dir = ps_cache.get_ligand_dir(scan_id, ligand_name)
    ligand_path = ligand_dir / f"{ligand_name}.sdf"
    if not ligand_path.exists():
        datastore.download_object(ligand_address, str(ligand_path))
    return str(ligand_path)


def _load_gene_pdbs_df(
    scan_id: str,
    gene_name: str,
) -> pd.DataFrame:
    """
    Load the per-gene PDB table produced by the PDB preparation step.

    Parameters
    ----------
    scan_id : str
        Identifier for the scan cache namespace.
    gene_name : str
        Gene symbol used to locate the per-gene cache folder.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by PDB id with metadata and cached file paths.

    Examples
    --------
    >>> df = _load_gene_pdbs_df("scan1", "TP53")
    """
    csv_path = ps_cache.gene_pdbs_csv_path(scan_id, gene_name)
    if not csv_path.exists():
        raise FileNotFoundError(f"pdbs file not found for {gene_name}: {csv_path}. "
                                "Run pdb_clean first.")
    return pd.read_csv(csv_path, index_col="id")


def _resolve_pdb_path(row_path: str, scan_id: str, gene_name: str) -> str:
    """
    Resolve a usable PDB path from a cached table row.

    Parameters
    ----------
    row_path : str
        Path value from the CSV row.
    scan_id : str
        Identifier for the scan cache namespace.
    gene_name : str
        Gene symbol used to locate the per-gene cache folder.

    Returns
    -------
    str
        Resolved path. This may still be non-existent if inputs are invalid.

    Examples
    --------
    >>> _resolve_pdb_path("cleaned.pdb", "scan1", "TP53")
    """
    if row_path and os.path.isabs(row_path) and os.path.exists(row_path):
        return row_path
    gene_dir = ps_cache.get_gene_dir(scan_id, gene_name)
    candidate = gene_dir / os.path.basename(row_path) if row_path else None
    if candidate and candidate.exists():
        return str(candidate)
    if row_path and os.path.exists(row_path):
        return row_path
    return str(candidate) if candidate else row_path


def _write_complex_pdb(complex_tuple, output_path: str) -> bool:
    """
    Write a docked complex to disk as a PDB file.

    Parameters
    ----------
    complex_tuple : object
        Complex object returned by the docking backend.
    output_path : str
        Destination PDB path.

    Returns
    -------
    bool
        True if the file was written successfully, otherwise False.

    Examples
    --------
    >>> ok = _write_complex_pdb(([], []), "out.pdb")
    """
    try:
        from rdkit import Chem  # type: ignore
    except ImportError:
        raise ImportError("RDKit is required for run_docking but not installed")
    try:
        complex_mol = Chem.CombineMols(complex_tuple[0][0], complex_tuple[0][1])
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        Chem.rdmolfiles.MolToPDBFile(complex_mol, output_path)
        return True
    except Exception as e:  # noqa: BLE001
        log_progress("run_docking", 0, f"Failed to write complex PDB: {e}")
        return False


def _upload_file(path: str, key: str, file_type: str, data_type: str) -> str:
    """
    Upload a local file to the datastore.

    Parameters
    ----------
    path : str
        Path to the local file.
    key : str
        Destination key or name in the datastore.
    file_type : str
        File type metadata for the DataCard.
    data_type : str
        Data type metadata for the DataCard.

    Returns
    -------
    str
        Datastore address of the uploaded object.

    Examples
    --------
    >>> addr = _upload_file("file.txt", "k", "txt", "text/plain")
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")
    with open(path, "rb") as f:
        content = f.read()
    card = DataCard(address="", file_type=file_type, data_type=data_type)
    if file_type in ("json", "csv", "pdb", "sdf"):
        # text upload
        with open(path, "r") as f:
            text_content = f.read()
        return datastore.upload_data_from_memory(text_content, key, card)
    return datastore.upload_data_from_memory(content, key, card)


def run_docking(
    gene_name: str,
    ligand_name: str,
    ligand_address: str,
    scan_id: str,
    output: str,
    exhaustiveness: int = 32,
    num_modes: int = 8,
) -> str:
    """
    Run docking for a single gene and ligand pair.

    Parameters
    ----------
    gene_name : str
        Gene symbol. The PDB preparation step must already have been run for
        this gene under the same scan_id.
    ligand_name : str
        Ligand label used for namespacing local outputs.
    ligand_address : str
        Datastore address of the prepared ligand SDF.
    scan_id : str
        Scan identifier used to locate cached inputs and group outputs.
    output : str
        Output name prefix for uploaded datastore artifacts.
    exhaustiveness : int, default 32
        Docking search exhaustiveness.
    num_modes : int, default 8
        Maximum number of poses to generate per target structure.

    Returns
    -------
    str
        DeepChem address of a summary JSON describing the uploaded
        artifacts and the local cache paths.

    Examples
    --------
    >>> addr = run_docking(
    ...     gene_name="TP53",
    ...     ligand_name="LIG",
    ...     ligand_address="ds://ligand.sdf",
    ...     scan_id="scan1",
    ...     output="my_run",
    ... )
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")
    if not gene_name or not ligand_name:
        raise ValueError("gene_name and ligand_name are required")
    if not ligand_address:
        raise ValueError("ligand_address is required")
    if not scan_id:
        raise ValueError("scan_id is required")
    if not output:
        raise ValueError("output is required")

    try:
        from rdkit import Chem  # type: ignore  # noqa: F401
    except ImportError:
        raise ImportError("RDKit is required for run_docking but not installed")

    top_csv_path = ps_cache.top_score_gene_ligand_csv_path(scan_id, gene_name, ligand_name)
    complexes_dir = ps_cache.get_gene_complexes_dir(scan_id, gene_name)

    if top_csv_path.exists():
        log_progress("run_docking", 10, f"Skipping {gene_name} {ligand_name}: {top_csv_path} exists")
        sorted_maindf = pd.read_csv(top_csv_path)
    else:
        log_progress("run_docking", 5, f"Preparing docking for {gene_name} {ligand_name}")
        df = _load_gene_pdbs_df(scan_id, gene_name)

        log_progress("run_docking", 10, "Downloading ligand to cache")
        ligand_path = _download_ligand_to_cache(ligand_address, scan_id, ligand_name)

        work_dir = tempfile.mkdtemp(prefix=f"dock_{gene_name}_")
        try:
            docking_scores: Dict[str, List[Optional[float]]] = {}
            for i, pdb_id in enumerate(list(df.index)):
                row = df.loc[pdb_id]
                pdb_path = _resolve_pdb_path(str(row["path"]), scan_id, gene_name)
                score: List[Optional[float]] = [None]
                complex = None
                try:
                    assert os.path.exists(pdb_path), f"pdb file not found for {pdb_id}: {pdb_path}"
                    assert os.path.exists(ligand_path), f"ligand file not found: {ligand_path}"
                    complex, score = _vina_docking(  # type: ignore
                        pdb_path,
                        ligand_path,
                        work_dir,
                        exhaustiveness=exhaustiveness,
                        num_modes=num_modes,
                    )
                except Exception as e:  # noqa: BLE001
                    log_progress("run_docking", 0, f"Docking failed for {pdb_id} => {e}")
                if complex is not None:
                    out_pdb = ps_cache.complex_path(scan_id, gene_name, str(pdb_id), ligand_name)
                    _write_complex_pdb(complex, str(out_pdb))
                docking_scores[str(pdb_id)] = list(score)
                progress = 20 + int(60 * (i + 1) / max(1, len(df)))
                log_progress("run_docking", progress, f"docked {pdb_id} ({i + 1}/{len(df)})")
        finally:
            import shutil

            shutil.rmtree(work_dir, ignore_errors=True)

        df_docked = pd.DataFrame({
            "id": list(docking_scores.keys()),
            "scores": list(docking_scores.values()),
        })
        df_docked = df_docked.dropna()
        df_docked["top_score"] = df_docked["scores"].apply(lambda x: x[0] if x else None)
        df_docked = df_docked.set_index("id")

        merged_df = pd.merge(df, df_docked, how="left", left_index=True, right_index=True)
        top_df = merged_df[["chains", "resolution", "coverage", "top_score", "scores"]].sort_values(by="top_score")
        top2_df = top_df.iloc[[0, 1]] if len(top_df) > 2 else top_df
        top2_df = top2_df.copy()
        top2_df["gene_name"] = [gene_name] * len(top2_df)

        sorted_maindf = top2_df.sort_values(by="top_score")
        os.makedirs(str(top_csv_path.parent), exist_ok=True)
        sorted_maindf.to_csv(top_csv_path)

    log_progress("run_docking", 90, "uploading docking artifacts to datastore")
    csv_card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    csv_key = f"{output}_top_score_{gene_name}_{ligand_name}.csv"
    df_to_upload = pd.read_csv(top_csv_path)
    top_csv_address = datastore.upload_data_from_memory(df_to_upload, csv_key, csv_card)

    complex_addresses: Dict[str, str] = {}
    pdb_card = DataCard(address="", file_type="pdb", data_type="text/plain")
    for complex_file in sorted(complexes_dir.glob(f"complex_{gene_name}_*_{ligand_name}.pdb")):
        with open(complex_file, "r") as f:
            pdb_content = f.read()
        key = f"{output}_{complex_file.stem}.pdb"
        complex_addresses[complex_file.stem] = datastore.upload_data_from_memory(pdb_content, key, pdb_card)

    summary = {
        "gene_name": gene_name,
        "ligand_name": ligand_name,
        "scan_id": scan_id,
        "top_score_csv_address": top_csv_address,
        "top_score_csv_path": str(top_csv_path),
        "complexes_dir": str(complexes_dir),
        "complex_addresses": complex_addresses,
        "exhaustiveness": exhaustiveness,
        "num_modes": num_modes,
    }
    summary_json = json.dumps(summary, indent=2, default=str)
    json_card = DataCard(address="", file_type="json", data_type="json")
    summary_key = f"{output}_run_docking_{gene_name}_{ligand_name}.json"
    summary_address = datastore.upload_data_from_memory(summary_json, summary_key, json_card)

    log_progress("run_docking", 100, "run_docking complete")
    return summary_address


__all__ = ["run_docking"]
