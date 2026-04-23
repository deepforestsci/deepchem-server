"""
Combine per-gene docking results into per-ligand summary tables.
"""

import json
import os
from typing import Dict, List

import pandas as pd

from deepchem_server.core.common import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.common.progress_logger import log_progress
from deepchem_server.core.primitives.proteome_scan import cache as ps_cache


def _concat_csv_from_folder(folder_path: str) -> pd.DataFrame:
    """
    Concatenate all CSV files in a folder into a single DataFrame.

    Parameters
    ----------
    folder_path : str
        Directory containing CSV files.

    Returns
    -------
    pandas.DataFrame
        Concatenated table. If the directory is missing or contains no CSV
        files, an empty DataFrame is returned.

    Examples
    --------
    >>> df = _concat_csv_from_folder("/tmp/does_not_exist")
    """
    main_df = pd.DataFrame()
    if not os.path.isdir(folder_path):
        return main_df
    for filename in sorted(os.listdir(folder_path)):
        file_path = os.path.join(folder_path, filename)
        if filename.endswith(".csv") and os.path.isfile(file_path):
            temp_df = pd.read_csv(file_path)
            main_df = pd.concat([main_df, temp_df], ignore_index=True)
    return main_df


def parse_results(
    scan_id: str,
    ligands: List[str],
    output: str,
) -> str:
    """Aggregate docking results and upload per-ligand tables.

    Parameters
    ----------
    scan_id: str
        Scan identifier. Results for each ligand are read from
        ``<cache_root>/<scan_id>/<ligand>/`` and written to
        ``<cache_root>/<scan_id>/scan_results/``.
    ligands: List[str]
        Ligand names to aggregate. Each must have a corresponding
        on-disk folder populated by ``run_docking``.
    output: str
        Output name prefix for uploaded datastore artifacts.

    Returns
    -------
    str
        DeepChem address of a summary JSON describing the aggregated
        per-ligand CSVs (address + local path).

    Examples
    --------
    >>> addr = parse_results("scan1", ["LIG"], output="my_run")
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")
    if not scan_id:
        raise ValueError("scan_id is required")
    if not ligands:
        raise ValueError("ligands is required")
    if not output:
        raise ValueError("output is required")

    ps_cache.get_scan_results_dir(scan_id)

    aggregated: Dict[str, Dict[str, object]] = {}
    n_ligands = len(ligands)
    for i, ligand in enumerate(ligands):
        raw_results_folder = str(ps_cache.get_ligand_dir(scan_id, ligand))
        log_progress(
            "parse_results", int(5 + 60 * i / max(1, n_ligands)), f"aggregating ligand {ligand}"
        )
        df = _concat_csv_from_folder(raw_results_folder)
        if df.empty:
            log_progress("parse_results", 0, f"no CSVs found for ligand {ligand}")
            aggregated[ligand] = {
                "address": None,
                "local_path": None,
                "num_rows": 0,
            }
            continue
        df = df.sort_values(by="top_score")
        if "gene_name" in df.columns:
            df = df.drop_duplicates("gene_name", keep="first")
        df = df.rename(columns={"Unnamed: 0": "pdb_id"})

        out_path = ps_cache.top_score_ligand_csv_path(scan_id, ligand)
        df.to_csv(out_path, index=False)

        csv_card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
        csv_key = f"{output}_top_score_{ligand}.csv"
        csv_address = datastore.upload_data_from_memory(df, csv_key, csv_card)
        aggregated[ligand] = {
            "address": csv_address,
            "local_path": str(out_path),
            "num_rows": int(len(df)),
        }

    summary = {
        "scan_id": scan_id,
        "ligands": ligands,
        "aggregated": aggregated,
    }
    summary_json = json.dumps(summary, indent=2, default=str)
    json_card = DataCard(address="", file_type="json", data_type="json")
    summary_address = datastore.upload_data_from_memory(
        summary_json, f"{output}_parse_results.json", json_card
    )

    log_progress("parse_results", 100, "parse_results complete")
    return summary_address


__all__ = ["parse_results"]
