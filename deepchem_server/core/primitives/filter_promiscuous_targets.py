import json
import os
from typing import Dict, List

import pandas as pd

from deepchem_server.core.common import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.common.progress_logger import log_progress


def filter_promiscuous_targets(
    scan_result_addresses: List[str],
    thresholds: List[List],
    output: str,
) -> str:
    """
    Filter the promiscuous targets for the given thresholds and scan results directory

    Parameters
    ----------
    scan_result_addresses : List[str]
        DeepChem addresses of scan result CSV files.
    thresholds : List[List]
        List of [m, n] pairs where m is the percentile cutoff (0-100) and
        n is the minimum occurrence count across ligands.
    output : str
        Output name prefix used for all uploaded result files.

    Returns
    -------
    str
        DeepChem address of the JSON results file containing:

        - promiscuous_targets: dictionary mapping threshold string to list of promiscuous gene names
        - promiscuous_targets_address: DeepChem address of the standalone promiscuous targets JSON file
        - filtered_results: dictionary mapping threshold string to a dictionary of {original_filename: filtered_csv_address}

    Raises
    ------
    ValueError
        If the datastore is not configured or no scan result addresses are provided.
    KeyError
        If any CSV is missing the required gene_name column.

    Examples
    --------
    >>> result_address = filter_promiscuous_targets(
    ...     scan_result_addresses=[
    ...         "deepchem://user/project/LigandA.csv",
    ...         "deepchem://user/project/LigandB.csv",
    ...     ],
    ...     thresholds=[[15, 1]],
    ...     output="promiscuity_run1",
    ... )
    >>> print(result_address)
    deepchem://user/project/promiscuity_run1_results.json
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")

    if not scan_result_addresses:
        raise ValueError("At least one scan result address is required")

    log_progress("filter_promiscuous_targets", 5, "downloading scan result CSVs from datastore")

    # Download all ligand DataFrames with their address filename as the key
    dfs: Dict[str, pd.DataFrame] = {}
    for address in scan_result_addresses:
        df = datastore.get(address)
        if not isinstance(df, pd.DataFrame):
            raise ValueError(f"Expected CSV at address {address}, got {type(df)}")
        if "gene_name" not in df.columns:
            raise KeyError(f"CSV at {address} is missing required 'gene_name' column")
        filename = os.path.basename(address.split("://")[-1])
        dfs[filename] = df

    log_progress("filter_promiscuous_targets", 20, "identifying promiscuous targets for each threshold")

    promiscuous_targets: Dict[str, List[str]] = {}
    for threshold_pair in thresholds:
        m, n = int(threshold_pair[0]), int(threshold_pair[1])
        target_threshold = m / 100.0

        top_slices: List[pd.DataFrame] = []
        for df in dfs.values():
            top_slices.append(df.iloc[:int(len(df) * target_threshold) + 1])

        gene_counts = (pd.concat(top_slices)["gene_name"].value_counts().rename_axis("gene_name").reset_index(
            name="occurrence in top"))
        promiscuous_genes: List[str] = gene_counts[gene_counts["occurrence in top"] >= n]["gene_name"].tolist()

        key = f"{m}%_{n}"
        promiscuous_targets[key] = promiscuous_genes

    log_progress("filter_promiscuous_targets", 40, "uploading promiscuous targets JSON")

    promiscuous_targets_json = json.dumps(promiscuous_targets, indent=2)
    promiscuous_card = DataCard(address="", file_type="json", data_type="json")
    promiscuous_targets_address = datastore.upload_data_from_memory(
        promiscuous_targets_json,
        f"{output}_promiscuous_targets.json",
        promiscuous_card,
    )

    log_progress("filter_promiscuous_targets", 50, "filtering and uploading per-threshold CSVs")

    filtered_results: Dict[str, Dict[str, str]] = {}
    n_thresholds = len(thresholds)
    for i, threshold_pair in enumerate(thresholds):
        m, n = int(threshold_pair[0]), int(threshold_pair[1])
        key = f"{m}%_{n}"
        genes_to_remove = promiscuous_targets[key]

        filtered_addresses: Dict[str, str] = {}
        for filename, df in dfs.items():
            filtered_df = df[~df["gene_name"].isin(genes_to_remove)]
            csv_name = f"{output}_filtered_{key}_{filename}"
            csv_card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
            filtered_address = datastore.upload_data_from_memory(filtered_df, csv_name, csv_card)
            filtered_addresses[filename] = filtered_address

        filtered_results[key] = filtered_addresses
        progress = 50 + int(40 * (i + 1) / n_thresholds)
        log_progress("filter_promiscuous_targets", progress, f"uploaded filtered CSVs for threshold {key}")

    log_progress("filter_promiscuous_targets", 90, "uploading results summary")

    results = {
        "promiscuous_targets": promiscuous_targets,
        "promiscuous_targets_address": promiscuous_targets_address,
        "filtered_results": filtered_results,
    }
    results_json = json.dumps(results)
    results_card = DataCard(address="", file_type="json", data_type="json")
    result_address = datastore.upload_data_from_memory(
        results_json,
        f"{output}_results.json",
        results_card,
    )

    log_progress("filter_promiscuous_targets", 100, "filter_promiscuous_targets completed successfully")
    return result_address
