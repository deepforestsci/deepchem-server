"""DEL denoising primitive.

Two scoring modes:

1. Unified - Poisson confidence-interval ratio across all replicates.
2. Non-unified - z-score computed independently for target and control.

Optionally collapses three-part rows into pairwise combinations before
scoring (use_disynthon_pairs=True).
"""

import logging
import os
import json
import tempfile
from math import sqrt
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2

from deepchem_server.core.common import config
from deepchem_server.core.common.address import DeepchemAddress
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.common.progress_logger import log_progress


logger = logging.getLogger(__name__)

DEFAULT_SMILES_COLS = ["smiles_a", "smiles_b", "smiles_c"]
DEFAULT_CONTROL_COLS = ["seq_matrix_1", "seq_matrix_2", "seq_matrix_3"]
DEFAULT_TARGET_COLS = ["seq_target_1", "seq_target_2", "seq_target_3"]


def _poissfit(vec: pd.Series, alpha: float = 0.05) -> Tuple[float, float]:
    """Poisson confidence interval for replicate counts.

    Parameters
    ----------
    vec : pd.Series
        Replicate counts for one row.
    alpha : float
        Significance level (default 0.05 for 95% CI).

    Returns
    -------
    Tuple[float, float]
        (lower_bound, upper_bound) of the estimated rate.
    """
    k_sum = vec.sum()
    n = len(vec)
    lower = 0.5 * chi2.ppf(alpha / 2, 2 * k_sum) / n
    upper = 0.5 * chi2.ppf(1 - alpha / 2, 2 * (k_sum + 1)) / n
    return (lower, upper)


def _get_enrichment_ratio(row: pd.Series,
                          control_cols: List[str],
                          target_cols: List[str],
                          alpha: float = 0.05) -> float:
    """Enrichment ratio: target_lower_bound / control_upper_bound.

    Parameters
    ----------
    row : pd.Series
        One row with control and target count columns.
    control_cols : List[str]
        Control count column names.
    target_cols : List[str]
        Target count column names.
    alpha : float
        Significance level for the confidence interval.

    Returns
    -------
    float
        Ratio or 0.0 if control upper bound is zero.
    """
    _, c_upper = _poissfit(row[control_cols], alpha)
    t_lower, _ = _poissfit(row[target_cols], alpha)
    if c_upper == 0:
        return 0.0
    return t_lower / c_upper


def _calculate_poisson_enrichment(df: pd.DataFrame,
                                  control_cols: List[str],
                                  target_cols: List[str],
                                  alpha: float = 0.05) -> pd.DataFrame:
    """Add a Poisson_Enrichment column to the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Input data with control and target count columns.
    control_cols : List[str]
        Control count column names.
    target_cols : List[str]
        Target count column names.
    alpha : float
        Significance level for confidence intervals.

    Returns
    -------
    pd.DataFrame
        Copy of dataframe with a Poisson_Enrichment column added.
    """
    result_df = df.copy()
    sub_df = result_df[control_cols + target_cols].astype(float)
    result_df["Poisson_Enrichment"] = sub_df.apply(
        lambda row: _get_enrichment_ratio(row, control_cols, target_cols, alpha), axis=1)
    return result_df


def _calculate_normalized_enrichment_score(row: pd.Series, total_sum: float, row_count: int, column_name: str) -> float:
    """Z-score for one row: (p0 - p1) / sqrt(p1 * (1 - p1)).

    Parameters
    ----------
    row : pd.Series
        One DataFrame row.
    total_sum : float
        Sum of column_name across all rows.
    row_count : int
        Number of rows in the DataFrame.
    column_name : str
        Column to read the count from.

    Returns
    -------
    float
        Normalized score.
    """
    p0 = row[column_name] / total_sum
    p1 = 1 / row_count
    return (p0 - p1) / sqrt(p1 * (1 - p1))


def _calculate_hit_threshold(df: pd.DataFrame, column_name: str, percentile: float) -> float:
    """Return the percentile cutoff for a column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    column_name : str
        Column to compute the percentile on.
    percentile : float
        Percentile value (0--100).

    Returns
    -------
    float
        The cutoff value.
    """
    return np.percentile(df[column_name], percentile)


def _get_disynthon_smiles(
    d1_idx: str,
    d2_idx: str,
    smiles_dict_inv: Dict[str, str],
    failed_smiles: Set,
    failed_combines: Set,
) -> Optional[str]:
    """Look up two fragments by index and merge them into one SMILES string.

    Parameters
    ----------
    d1_idx : str
        Index key for the first fragment.
    d2_idx : str
        Index key for the second fragment.
    smiles_dict_inv : Dict[str, str]
        Maps index keys to SMILES strings.
    failed_smiles : Set
        Collects invalid SMILES (modified in place).
    failed_combines : Set
        Collects pairs that fail to merge (modified in place).

    Returns
    -------
    Optional[str]
        Merged SMILES or None on failure.
    """
    from rdkit import Chem

    try:
        smi_1 = smiles_dict_inv[d1_idx]
        smi_2 = smiles_dict_inv[d2_idx]
    except KeyError:
        return None

    mol1 = Chem.MolFromSmiles(smi_1)
    if mol1 is None:
        failed_smiles.add(smi_1)
        return None
    mol2 = Chem.MolFromSmiles(smi_2)
    if mol2 is None:
        failed_smiles.add(smi_2)
        return None
    try:
        combined = Chem.CombineMols(mol1, mol2)
        return Chem.MolToSmiles(combined)
    except Exception:
        failed_combines.add((smi_1, smi_2))
        return None


def _create_disynthon_pairs(
    df: pd.DataFrame,
    smiles_cols: List[str],
    count_cols: List[str],
    is_unified: bool,
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """Generate all pairwise groupings (AB, AC, BC) from three-part data.

    Parameters
    ----------
    df : pd.DataFrame
        Input data with SMILES and count columns.
    smiles_cols : List[str]
        Three SMILES column names.
    count_cols : List[str]
        Count columns to aggregate.
    is_unified : bool
        If True, keep individual count columns.  If False,
        pre-sum them into two totals before grouping.

    Returns
    -------
    Tuple[pd.DataFrame, Dict[str, str]]
        (pair_df, smiles_dict).  pair_df has Disynthon_1,
        Disynthon_2 and aggregated counts.  smiles_dict maps
        SMILES to index strings.
    """
    smiles_set: set = set()
    for col in smiles_cols:
        smiles_set.update(df[col].dropna())
    smiles_list = list(smiles_set)
    smiles_dict = {smi: str(i) for i, smi in enumerate(smiles_list)}

    df_work = df.copy()
    for col in smiles_cols:
        df_work[col] = df_work[col].map(smiles_dict)
    df_work = df_work.dropna(subset=smiles_cols)

    if not is_unified:
        target_count_cols = [c for c in count_cols if "target" in c.lower()]
        control_count_cols = [c for c in count_cols if "matrix" in c.lower() or "control" in c.lower()]
        df_work["seq_target_sum"] = df_work[target_count_cols].sum(axis=1)
        df_work["seq_control_sum"] = df_work[control_count_cols].sum(axis=1)
        agg_cols = ["seq_target_sum", "seq_control_sum"]
    else:
        agg_cols = count_cols

    pair_frames = []
    n = len(smiles_cols)
    for i in range(n):
        for j in range(i + 1, n):
            col1, col2 = smiles_cols[i], smiles_cols[j]
            agg_dict = {c: "sum" for c in agg_cols}
            pair_df = df_work.groupby([col1, col2]).agg(agg_dict).reset_index()
            pair_df.rename(columns={col1: "Disynthon_1", col2: "Disynthon_2"}, inplace=True)
            pair_frames.append(pair_df)

    result = pd.concat(pair_frames, ignore_index=True)
    return result, smiles_dict


def _collapse_to_disynthons(
    df: pd.DataFrame,
    smiles_cols: List[str],
    control_cols: List[str],
    target_cols: List[str],
    is_unified: bool,
    aggregate_operation: str = "sum",
    min_count_threshold: int = 0,
) -> Tuple[pd.DataFrame, int]:
    """Collapse three-part rows into pairwise combinations.

    Parameters
    ----------
    df : pd.DataFrame
        Cleaned input (no NaN/duplicate rows).
    smiles_cols : List[str]
        Three SMILES column names.
    control_cols : List[str]
        Control count column names.
    target_cols : List[str]
        Target count column names.
    is_unified : bool
        If True, keep individual count columns.  If False,
        pre-sum into two totals.
    aggregate_operation : str
        How to combine duplicate counts: 'sum' or 'mean'.
    min_count_threshold : int
        Drop rows with total count below this value.

    Returns
    -------
    Tuple[pd.DataFrame, int]
        (collapsed_df, n_failed).  collapsed_df has a disynthons
        column and aggregated counts.  n_failed is the number of
        SMILES that could not be merged.
    """
    count_cols = control_cols + target_cols

    pair_df, smiles_dict = _create_disynthon_pairs(df, smiles_cols, count_cols, is_unified)

    smiles_dict_inv = {v: k for k, v in smiles_dict.items()}
    failed_smiles: Set = set()
    failed_combines: Set = set()

    pair_df["disynthons"] = pair_df.apply(
        lambda row: _get_disynthon_smiles(row["Disynthon_1"], row["Disynthon_2"], smiles_dict_inv, failed_smiles,
                                          failed_combines),
        axis=1,
    )

    pair_df = pair_df[pair_df["disynthons"].notna()]
    n_failed = len(failed_smiles) + len(failed_combines)
    logger.info(f"del_denoise: disynthon collapse had {n_failed} SMILES combination failures")

    numeric_cols = [c for c in pair_df.columns if c not in ("Disynthon_1", "Disynthon_2", "disynthons")]
    agg_dict = {c: aggregate_operation for c in numeric_cols}
    pair_df = pair_df.groupby("disynthons").agg(agg_dict).reset_index()

    if min_count_threshold > 0:
        count_sum = pair_df[numeric_cols].sum(axis=1)
        pair_df = pair_df[count_sum >= min_count_threshold]

    return pair_df, n_failed


def del_denoise(
    dataset_address: str,
    output_key: str,
    strategy: str = "unified",
    control_cols: Optional[List[str]] = None,
    target_cols: Optional[List[str]] = None,
    add_hit_labels: bool = False,
    hit_percentile: float = 90.0,
    alpha: float = 0.05,
    drop_duplicates: bool = True,
    use_disynthon_pairs: bool = False,
    smiles_cols: Optional[List[str]] = None,
    aggregate_operation: str = "sum",
    min_count_threshold: int = 0,
) -> str:
    """Score DEL screening data to find strong binders.

    Reads a CSV of raw counts, scores each row and writes the
    result back to the datastore.  Two strategies are available:
    'unified' (Poisson-based) and 'non_unified' (z-score-based).

    Parameters
    ----------
    dataset_address : str
        Datastore address of the input CSV.
    output_key : str
        Name for the output CSV in the datastore.
    strategy : str
        'unified' (Poisson ratio) or 'non_unified' (z-score).
    control_cols : Optional[List[str]]
        Control count column names.
    target_cols : Optional[List[str]]
        Target count column names.
    add_hit_labels : bool
        Add binary 0/1 hit columns based on a percentile cutoff.
    hit_percentile : float
        Percentile cutoff for hits (0--100).  Used when
        add_hit_labels is True.
    alpha : float
        Confidence level for Poisson intervals.  Used when
        strategy is 'unified'.
    drop_duplicates : bool
        Remove duplicate SMILES rows before scoring.
    use_disynthon_pairs : bool
        Collapse three-part rows into pairwise combinations
        before scoring.
    smiles_cols : Optional[List[str]]
        Three SMILES column names for the pairwise collapse.
        Used when use_disynthon_pairs is True.
    aggregate_operation : str
        'sum' or 'mean' for combining duplicate pair counts.
        Used when use_disynthon_pairs is True.
    min_count_threshold : int
        Drop pair rows with total count below this value.
        Used when use_disynthon_pairs is True.

    Returns
    -------
    str
        Datastore address of the output CSV.

    Raises
    ------
    ValueError
        If strategy is invalid or the datastore is not configured.

    Examples
    --------
    Unified scoring:

    >>> from deepchem_server.core.common.cards import DataCard
    >>> from deepchem_server.core.common import config
    >>> from deepchem_server.core.datastore import DiskDataStore
    >>> import tempfile, pandas as pd
    >>> disk_datastore = DiskDataStore('profile', 'project', tempfile.mkdtemp())
    >>> config.set_datastore(disk_datastore)
    >>> df = pd.DataFrame({
    ...     "smiles": ["CCO", "CCN", "CCC"],
    ...     "seq_matrix_1": [10, 20, 5], "seq_matrix_2": [12, 18, 6],
    ...     "seq_matrix_3": [11, 22, 4],
    ...     "seq_target_1": [50, 30, 8], "seq_target_2": [55, 28, 7],
    ...     "seq_target_3": [48, 32, 9],
    ... })
    >>> card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    >>> addr = disk_datastore.upload_data_from_memory(df, "raw_del.csv", card)
    >>> result_addr = del_denoise(dataset_address=addr, output_key="denoised")
    >>> result_addr
    'deepchem://project/profile/denoised.csv'

    With hit labels:

    >>> result_addr = del_denoise(
    ...     dataset_address=addr,
    ...     output_key="denoised_hits",
    ...     strategy="unified",
    ...     add_hit_labels=True,
    ...     hit_percentile=90.0,
    ... )
    >>> result_addr
    'deepchem://project/profile/denoised_hits.csv'

    Non-unified scoring:

    >>> result_addr = del_denoise(
    ...     dataset_address=addr,
    ...     output_key="denoised_nu",
    ...     strategy="non_unified",
    ...     add_hit_labels=True,
    ... )
    >>> result_addr
    'deepchem://project/profile/denoised_nu.csv'
    """
    if control_cols is None:
        control_cols = list(DEFAULT_CONTROL_COLS)
    if target_cols is None:
        target_cols = list(DEFAULT_TARGET_COLS)
    if smiles_cols is None:
        smiles_cols = list(DEFAULT_SMILES_COLS)

    # All params arrive as plain strings when called through the HTTP router.
    add_hit_labels = str(add_hit_labels).lower() == "true"
    drop_duplicates = str(drop_duplicates).lower() == "true"
    use_disynthon_pairs = str(use_disynthon_pairs).lower() == "true"
    hit_percentile = float(hit_percentile)
    alpha = float(alpha)
    min_count_threshold = int(min_count_threshold)
    if isinstance(control_cols, str):
        control_cols = json.loads(control_cols)
    if isinstance(target_cols, str):
        target_cols = json.loads(target_cols)
    if isinstance(smiles_cols, str):
        smiles_cols = json.loads(smiles_cols)

    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")

    log_progress("del_denoise", 10, "downloading dataset")
    tmpdir = tempfile.TemporaryDirectory()
    local_path = os.path.join(tmpdir.name, "input.csv")
    datastore.download_object(dataset_address, local_path)

    log_progress("del_denoise", 15, "loading and cleaning data")
    df = pd.read_csv(local_path)
    n_input = len(df)

    if use_disynthon_pairs:
        df = df.dropna(subset=["smiles"])
        if drop_duplicates:
            df = df.drop_duplicates(subset=["smiles"])

        is_unified = strategy == "unified"
        log_progress("del_denoise", 25, "collapsing trisynthons into disynthon pairs")
        df, n_failed = _collapse_to_disynthons(
            df,
            smiles_cols,
            control_cols,
            target_cols,
            is_unified,
            aggregate_operation,
            min_count_threshold,
        )

        if not is_unified:
            control_cols = ["seq_control_sum"]
            target_cols = ["seq_target_sum"]

        smiles_col = "disynthons"
    else:
        smiles_col = "disynthons" if "disynthons" in df.columns else "smiles"
        df = df.dropna(subset=[smiles_col])
        if drop_duplicates:
            df = df.drop_duplicates(subset=[smiles_col])

    n_output = len(df)
    logger.info(f"del_denoise: {n_input} -> {n_output} rows after cleaning")

    if strategy == "unified":
        log_progress("del_denoise", 50, "computing Poisson enrichment scores")
        df = _calculate_poisson_enrichment(df, control_cols, target_cols, alpha)

        if add_hit_labels:
            log_progress("del_denoise", 70, "computing hit labels")
            threshold = _calculate_hit_threshold(df, "Poisson_Enrichment", hit_percentile)
            df["hits"] = (df["Poisson_Enrichment"] > threshold).astype(int)

    elif strategy == "non_unified":
        log_progress("del_denoise", 40, "summing replicate counts")
        df["seq_target_sum"] = df[target_cols].sum(axis=1)
        df["seq_control_sum"] = df[control_cols].sum(axis=1)

        log_progress("del_denoise", 55, "computing z-score enrichment (target)")
        total_target = df["seq_target_sum"].sum()
        total_control = df["seq_control_sum"].sum()
        row_count = len(df)

        df["Target_Enrichment_Score"] = df.apply(
            lambda row: _calculate_normalized_enrichment_score(row, total_target, row_count, "seq_target_sum"),
            axis=1,
        )

        log_progress("del_denoise", 65, "computing z-score enrichment (control)")
        df["Control_Enrichment_Score"] = df.apply(
            lambda row: _calculate_normalized_enrichment_score(row, total_control, row_count, "seq_control_sum"),
            axis=1,
        )

        if add_hit_labels:
            log_progress("del_denoise", 75, "computing hit labels")
            target_threshold = _calculate_hit_threshold(df, "Target_Enrichment_Score", hit_percentile)
            control_threshold = _calculate_hit_threshold(df, "Control_Enrichment_Score", hit_percentile)
            df["target_hits"] = (df["Target_Enrichment_Score"] > target_threshold).astype(int)
            df["control_hits"] = (df["Control_Enrichment_Score"] > control_threshold).astype(int)
    else:
        raise ValueError(f"Unknown strategy '{strategy}'. Must be 'unified' or 'non_unified'.")

    log_progress("del_denoise", 90, "uploading denoised dataset")
    if not output_key.endswith(".csv"):
        output_key = output_key + ".csv"
    output_key = DeepchemAddress.get_key(output_key)

    card_kwargs: Dict = dict(
        address="",
        file_type="csv",
        data_type="pandas.DataFrame",
        shape=(n_output, len(df.columns)),
        description=f"DEL denoised enrichment scores ({strategy})",
        strategy=strategy,
        parent=dataset_address,
        control_cols=control_cols,
        target_cols=target_cols,
        n_input_rows=n_input,
        n_output_rows=n_output,
        use_disynthon_pairs=use_disynthon_pairs,
        add_hit_labels=str(add_hit_labels),
        hit_percentile=hit_percentile,
    )
    if use_disynthon_pairs:
        card_kwargs["n_disynthons"] = n_output
        card_kwargs["aggregate_operation"] = aggregate_operation
        card_kwargs["min_count_threshold"] = min_count_threshold
        card_kwargs["n_failed_smiles"] = n_failed

    if strategy == "unified":
        card_kwargs["alpha"] = alpha

    card = DataCard(**card_kwargs)

    return datastore.upload_data_from_memory(df, output_key, card)
