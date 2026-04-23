"""
Analyze docking poses and pocket features for the proteome scan workflow.
"""

import json
import multiprocessing
import os
import re
import shutil
import subprocess
from typing import Dict, List, Optional, Tuple

import pandas as pd

from deepchem_server.core.common import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.common.progress_logger import log_progress
from deepchem_server.core.primitives.proteome_scan import cache as ps_cache


def _import_pymol():
    """
    Import PyMOL modules used by this analysis.

    Returns
    -------
    tuple
        Tuple of (cmd, stored) objects from PyMOL.

    Examples
    --------
    >>> cmd, stored = _import_pymol()
    """
    try:
        from pymol import cmd, stored  # type: ignore
    except ImportError as e:
        raise ImportError("PyMOL is required for run_multi_pose_analysis but not installed") from e
    return cmd, stored


_FPOCKET_DOCKER_IMAGE = os.getenv("FPOCKET_DOCKER_IMAGE", "fpocket/fpocket")


def _run_fpocket(protein_path: str) -> None:
    """
    Run fpocket on a protein PDB using the Docker image.

    Parameters
    ----------
    protein_path : str
        Path to the protein PDB file. The parent directory is used as a shared
        working directory for the container outputs.

    Returns
    -------
    None

    Raises
    ------
    RuntimeError
        If the Docker container exits with a non-zero status.

    Examples
    --------
    >>> _run_fpocket("protein.pdb")
    """
    protein_path = os.path.abspath(protein_path)
    host_dir = os.path.dirname(protein_path)
    filename = os.path.basename(protein_path)

    try:
        os.chmod(host_dir, 0o777)
    except OSError:
        pass

    proc = subprocess.run(
        [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{host_dir}:/workdir",
            _FPOCKET_DOCKER_IMAGE,
            "fpocket",
            "-f",
            f"/workdir/{filename}",
        ],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"fpocket Docker container failed on {protein_path} (exit {proc.returncode}):\n{proc.stderr.strip()}"
        )


def _separate_protein_ligand(pose_path: str, run_dir: str) -> Tuple[str, str]:
    """
    Split a complex PDB into protein and ligand PDB files.

    Parameters
    ----------
    pose_path : str
        Path to the complex PDB file.
    run_dir : str
        Directory where the split files will be written.

    Returns
    -------
    tuple[str, str]
        Paths to (ligand_pdb_path, protein_pdb_path).

    Examples
    --------
    >>> lig, prot = _separate_protein_ligand("complex.pdb", "/tmp/run")
    """
    cmd, _ = _import_pymol()
    cmd.load(pose_path, "complex")
    cmd.select("ligand", "resn LIG+UNL")
    cmd.select("protein_1", "polymer and not resn LIG+UNL and not resn HOH")
    ligand_path = os.path.join(run_dir, "ligand.pdb")
    protein_path = os.path.join(run_dir, "protein.pdb")
    cmd.save(ligand_path, "ligand")
    cmd.save(protein_path, "protein_1")
    cmd.reinitialize()
    return ligand_path, protein_path


def _parse_pocket_data(file_path: str) -> list:
    """
    Parse fpocket pocket metadata from a protein_info.txt file.

    Parameters
    ----------
    file_path : str
        Path to the fpocket protein_info.txt file.

    Returns
    -------
    list
        List of pocket dictionaries.

    Examples
    --------
    >>> pockets = _parse_pocket_data("protein_info.txt")
    """
    with open(file_path, "r") as file:
        content = file.read()
    pocket_blocks = re.split(r"Pocket\s+\d+\s*:", content)
    pocket_headers = re.findall(r"Pocket\s+\d+\s*:", content)
    pockets = []
    for i, (_header, block) in enumerate(zip(pocket_headers, pocket_blocks[1:])):
        pocket_data: Dict[str, object] = {"pocket_id": i + 1}
        lines = block.strip().split("\n")
        for line in lines:
            if ":" in line:
                key, value = line.split(":", 1)
                key = key.strip()
                value = value.strip()
                try:
                    if "." in value:
                        value = float(value)
                    else:
                        value = int(value)
                except ValueError:
                    pass
                pocket_data[key] = value
        pockets.append(pocket_data)
    return pockets


def _analyse_overlaps(
    ligand_path: str,
    pockets_out_path: str,
    pockets_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute ligand and pocket overlap metrics using PyMOL.

    Parameters
    ----------
    ligand_path : str
        Path to the ligand PDB.
    pockets_out_path : str
        Path to the fpocket output PDB containing alpha spheres.
    pockets_df : pandas.DataFrame
        Pocket metadata table to augment with overlap metrics.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Full pockets table and the analysis table used for downstream metrics.

    Examples
    --------
    >>> full_df, analysis_df = _analyse_overlaps("ligand.pdb", "protein_out.pdb", pd.DataFrame())
    """
    cmd, stored = _import_pymol()
    cmd.load(ligand_path)
    cmd.load(pockets_out_path)
    stored.list = []
    cmd.iterate("(resn STP)", "stored.list.append(resi)")
    if not stored.list:
        cmd.reinitialize()
        raise ValueError("No STP residues found in fpocket output")
    last_stp = int(stored.list[-1])

    cmd.hide("lines", "resn STP")
    for my_index in range(1, last_stp + 1):
        pocket_name = f"pocket{my_index}"
        cmd.select(pocket_name, f"resn STP and resi {my_index}")
        cmd.show("spheres", pocket_name)
        cmd.set("sphere_scale", 0.3, pocket_name)
        cmd.set("sphere_transparency", 0.1, pocket_name)
        cmd.color(my_index, pocket_name)

    overlap_within_alpha_spheres: List[int] = []
    ratio_pocket_filled: List[float] = []
    for my_index in range(1, last_stp + 1):
        cmd.do(f"select overlapping{my_index}, ligand within 3 of pocket{my_index}")
        cmd.do(f"select p_overlapping{my_index}, pocket{my_index} within 3 of ligand")
        overlap_within_alpha_spheres.append(cmd.count_atoms(f"overlapping{my_index}"))
        pocket_alphasphere_count = cmd.count_atoms(f"pocket{my_index}")
        if pocket_alphasphere_count == 0:
            ratio_pocket_filled.append(0.0)
        else:
            ratio_pocket_filled.append(cmd.count_atoms(f"p_overlapping{my_index}") / pocket_alphasphere_count)

    VDWoverlap_data: List[float] = []
    for my_index in range(1, last_stp + 1):
        pocket_name = f"pocket{my_index}"
        VDWoverlap_data.append(cmd.overlap("ligand", pocket_name))

    ligand_atoms_count = cmd.count_atoms("ligand")
    if ligand_atoms_count == 0:
        cmd.reinitialize()
        raise ValueError("ligand has 0 atoms after loading")
    pockets_df["Ligand atoms overlap count"] = pd.Series(overlap_within_alpha_spheres)
    pockets_df["% Ligand inside pocket"] = pockets_df["Ligand atoms overlap count"] / ligand_atoms_count * 100
    pockets_df["% Pocket Filled"] = pd.Series(ratio_pocket_filled) * 100
    pockets_df["Ligand_VDWoverlap"] = pd.Series(VDWoverlap_data)
    analysis_df = pockets_df[
        [
            "pocket_id",
            "Ligand_VDWoverlap",
            "% Ligand inside pocket",
            "% Pocket Filled",
            "Ligand atoms overlap count",
            "Number of Alpha Spheres",
            "Druggability Score",
            "Pocket_Druggability_Rank",
            "Volume",
            "Mean local hydrophobic density",
        ]
    ]
    cmd.reinitialize()
    return pockets_df, analysis_df


def analyse_pose(pose_path: str, results_dir: str, is_clean_up: bool = False) -> None:
    """
    Run pose analysis for one complex PDB.

    Parameters
    ----------
    pose_path : str
        Path to the complex PDB file.
    results_dir : str
        Directory where per-complex CSV outputs will be written.
    is_clean_up : bool, default False
        If True, delete the per-complex working directory after completion.

    Returns
    -------
    None

    Examples
    --------
    >>> analyse_pose("complex.pdb", "/tmp/results", is_clean_up=True)
    """
    run_name = os.path.basename(pose_path).split(".")[0]
    run_dir = os.path.join(os.getcwd(), run_name)
    os.makedirs(run_dir, exist_ok=True)
    ligand_path, protein_path = _separate_protein_ligand(pose_path, run_dir)
    _run_fpocket(protein_path)
    pockets_out_path = os.path.join(run_dir, "protein_out", "protein_out.pdb")
    pockets_info_path = os.path.join(run_dir, "protein_out", "protein_info.txt")
    parsed_pockets = _parse_pocket_data(pockets_info_path)
    if not parsed_pockets:
        if is_clean_up:
            shutil.rmtree(run_dir, ignore_errors=True)
        raise ValueError("No pockets found!")
    pockets_df = pd.DataFrame(parsed_pockets)
    pockets_df["Pocket_Druggability_Rank"] = pockets_df["Druggability Score"].rank(ascending=False, method="min")
    pockets_df, analysis_df = _analyse_overlaps(ligand_path, pockets_out_path, pockets_df)
    final_df = analysis_df[analysis_df["Ligand_VDWoverlap"] > 0].sort_values(by=["Ligand_VDWoverlap"])
    output_dir = os.path.join(results_dir, run_name)
    os.makedirs(output_dir, exist_ok=True)
    pockets_df.to_csv(os.path.join(output_dir, "full_analysis.csv"), index=False)
    final_df.to_csv(os.path.join(output_dir, "result.csv"), index=False)
    if is_clean_up:
        shutil.rmtree(run_dir, ignore_errors=True)


def _run_analyse_pose_worker(args: Tuple[str, str, bool]) -> bool:
    """
    Worker wrapper for analyse_pose that returns success as a boolean.

    Parameters
    ----------
    args : tuple[str, str, bool]
        Tuple of (complex_path, results_dir, is_clean_up).

    Returns
    -------
    bool
        True if analysis completed successfully, otherwise False.

    Examples
    --------
    >>> ok = _run_analyse_pose_worker(("complex.pdb", "/tmp/results", True))
    """
    complex_path, results_dir, is_clean_up = args
    try:
        analyse_pose(pose_path=complex_path, results_dir=results_dir, is_clean_up=is_clean_up)
        return True
    except Exception as e:  # noqa: BLE001
        log_progress("run_multi_pose_analysis", 0, f"error {complex_path}: {e}")
        return False


def _get_overall_ligand_interactions(df: pd.DataFrame) -> float:
    """
    Compute a capped overall ligand interaction score from a per-pocket table.

    Parameters
    ----------
    df : pandas.DataFrame
        Pocket-level table with a percent ligand inside pocket column.

    Returns
    -------
    float
        Sum of per-pocket percentages, capped at 100.

    Examples
    --------
    >>> _get_overall_ligand_interactions(pd.DataFrame({"% Ligand inside pocket": [10, 20]}))
    30
    """
    percents = df["% Ligand inside pocket"]
    total = sum(percents)
    return total if total <= 100 else 100


def _get_total_top_n_bucket_percentages(df: pd.DataFrame, n: int) -> float:
    """
    Sum percent ligand inside pocket for the top-n ranked pockets.

    Parameters
    ----------
    df : pandas.DataFrame
        Pocket-level table with rank and percentage columns.
    n : int
        Rank cutoff. Pockets with rank less than or equal to n are included.

    Returns
    -------
    float
        Total percentage, capped at 100.

    Examples
    --------
    >>> df = pd.DataFrame({"Pocket_Druggability_Rank": [1, 2], "% Ligand inside pocket": [30.0, 40.0]})
    >>> _get_total_top_n_bucket_percentages(df, n=1)
    30.0
    """
    value = df[df["Pocket_Druggability_Rank"] <= n]["% Ligand inside pocket"].to_numpy().sum()
    return float(min(value, 100))


def _download_complex_to_cache(
    address: str,
    scan_id: Optional[str],
    output: str,
) -> str:
    """
    Download a complex PDB into a local cache directory.

    Parameters
    ----------
    address : str
        Datastore address for the complex PDB.
    scan_id : str or None
        If provided, download into the scan cache namespace.
    output : str
        Output name used to namespace downloads when scan_id is not provided.

    Returns
    -------
    str
        Local path to the downloaded complex PDB.

    Examples
    --------
    >>> path = _download_complex_to_cache("ds://complex.pdb", scan_id=None, output="run1")
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")
    if scan_id:
        base = ps_cache.get_scan_dir(scan_id) / "pose_analysis_inputs"
    else:
        base = ps_cache.get_cache_root() / f"pose_analysis_{output}"
    base.mkdir(parents=True, exist_ok=True)
    name = os.path.basename(address.split("://")[-1])
    if not name.endswith(".pdb"):
        name = f"{name}.pdb"
    local_path = base / name
    if not local_path.exists():
        datastore.download_object(address, str(local_path))
    return str(local_path)


def run_multi_pose_analysis(
    complex_addresses: List[str],
    output: str,
    num_processes: int = 4,
    is_clean_up: bool = True,
    scan_id: Optional[str] = None,
) -> str:
    """
    Run pose analysis across multiple docked complexes.

    Parameters
    ----------
    complex_addresses : list[str]
        Datastore addresses of complex PDB files, typically produced by docking.
    output : str
        Output name prefix for uploaded datastore artifacts.
    num_processes : int, default 4
        Number of parallel workers.
    is_clean_up : bool, default True
        If True, delete the per-complex working directory after each worker completes.
    scan_id : str or None, default None
        If provided, intermediates are downloaded under the scan cache namespace.

    Returns
    -------
    str
        DeepChem address of a summary JSON referencing the aggregated
        CSV and per-complex metrics.

    Examples
    --------
    >>> addr = run_multi_pose_analysis(
    ...     complex_addresses=["ds://complex_TP53_1ABC_LIG.pdb"],
    ...     output="my_run",
    ...     num_processes=1,
    ...     is_clean_up=True,
    ...     scan_id="scan1",
    ... )
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")
    if not complex_addresses:
        raise ValueError("complex_addresses is required")
    if not output:
        raise ValueError("output is required")

    log_progress("run_multi_pose_analysis", 5, "downloading complexes into cache")
    local_paths: List[str] = []
    for addr in complex_addresses:
        local_paths.append(_download_complex_to_cache(addr, scan_id, output))

    if scan_id:
        results_dir = str(ps_cache.get_scan_dir(scan_id) / f"pose_analysis_{output}")
    else:
        results_dir = str(ps_cache.get_cache_root() / f"pose_analysis_{output}" / "results")
    os.makedirs(results_dir, exist_ok=True)

    log_progress(
        "run_multi_pose_analysis",
        20,
        f"running analyse_pose on {len(local_paths)} complexes with {num_processes} workers",
    )

    worker_args = [(p, results_dir, is_clean_up) for p in local_paths]
    if num_processes <= 1:
        results = [_run_analyse_pose_worker(a) for a in worker_args]
    else:
        with multiprocessing.Pool(processes=num_processes) as pool:
            results = pool.map(_run_analyse_pose_worker, worker_args)

    failed_poses = [complex_addresses[i] for i, ok in enumerate(results) if not ok]

    data: List[Dict[str, object]] = []
    for i, complex_path in enumerate(local_paths):
        datapoint: Dict[str, object] = {}
        run_name = os.path.basename(complex_path).split(".")[0]
        datapoint["complex"] = run_name
        datapoint["address"] = complex_addresses[i]
        if not results[i]:
            datapoint["total % Ligand inside pockets"] = None
            for n in [1, 5, 10]:
                datapoint[f"total % Ligand inside pockets (top{n} pockets)"] = None
            data.append(datapoint)
            continue
        results_path = os.path.join(results_dir, run_name, "full_analysis.csv")
        if not os.path.exists(results_path):
            datapoint["total % Ligand inside pockets"] = None
            for n in [1, 5, 10]:
                datapoint[f"total % Ligand inside pockets (top{n} pockets)"] = None
            data.append(datapoint)
            continue
        df = pd.read_csv(results_path)
        datapoint["total % Ligand inside pockets"] = _get_overall_ligand_interactions(df)
        for n in [1, 5, 10]:
            datapoint[f"total % Ligand inside pockets (top{n} pockets)"] = _get_total_top_n_bucket_percentages(df, n)
        data.append(datapoint)

    df_results = pd.DataFrame(data)
    summary_csv_path = os.path.join(results_dir, "pose_analysis.csv")
    df_results.to_csv(summary_csv_path, index=False)

    log_progress("run_multi_pose_analysis", 85, "uploading aggregated CSV to datastore")
    csv_card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    csv_key = f"{output}_pose_analysis.csv"
    csv_address = datastore.upload_data_from_memory(df_results, csv_key, csv_card)

    summary = {
        "output": output,
        "scan_id": scan_id,
        "num_complexes": len(complex_addresses),
        "num_failed": len(failed_poses),
        "failed_addresses": failed_poses,
        "results_csv_address": csv_address,
        "results_csv_path": summary_csv_path,
        "results_dir": results_dir,
        "per_complex": data,
    }
    summary_json = json.dumps(summary, indent=2, default=str)
    json_card = DataCard(address="", file_type="json", data_type="json")
    summary_address = datastore.upload_data_from_memory(
        summary_json, f"{output}_run_multi_pose_analysis.json", json_card
    )

    log_progress("run_multi_pose_analysis", 100, "run_multi_pose_analysis complete")
    return summary_address


__all__ = [
    "run_multi_pose_analysis",
    "analyse_pose",
]
