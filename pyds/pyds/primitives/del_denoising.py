"""
DEL Denoise primitive module for DeepChem Server.

Contains the DelDenoise class for scoring DNA-encoded library screening data.
"""

from typing import Any, Dict, List, Optional

from .base import Primitive


class DelDenoise(Primitive):
    """
    Primitive for DEL denoising and enrichment scoring.

    Submits a denoising job to the DeepChem Server API's
    /primitive/del_denoise endpoint.  Two scoring strategies are
    available:

    * unified — Poisson confidence-interval ratio across all replicates.
    * non_unified — z-score computed independently for target and control.

    Example Usage
    -------------
    .. code-block:: python

       del_denoise = DelDenoise(settings)
       response = del_denoise.run(
           dataset_address="deepchem://project/profile/raw_del.csv",
           output_key="denoised",
           strategy="unified",
           add_hit_labels=True,
           hit_percentile=90.0,
       )
    """

    def run(
        self,
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
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Run the DEL denoise primitive.

        Args:
            dataset_address: Datastore address of the input CSV with raw counts.
            output_key: Name for the output CSV in the datastore.
            strategy: 'unified' (Poisson ratio) or 'non_unified'
                (z-score).  Defaults to 'unified'.
            control_cols: Control replicate count column names.
            target_cols: Target replicate count column names.
            add_hit_labels: Append a binary hit column based on a percentile
                cutoff of the enrichment score.
            hit_percentile: Percentile cutoff for hit labelling (0--100).
            alpha: Significance level for Poisson CIs (unified only).
            drop_duplicates: Remove duplicate SMILES rows before scoring.
            use_disynthon_pairs: Collapse trisynthon rows into pairwise
                disynthon combinations before scoring.
            smiles_cols: Three SMILES column names for the pairwise collapse.
            aggregate_operation: 'sum' or 'mean' for combining
                duplicate pair counts (disynthon mode only).
            min_count_threshold: Drop disynthon rows below this total count.
            profile_name: Profile name (uses settings if not provided).
            project_name: Project name (uses settings if not provided).

        Returns:
            Response dict containing the scored output address as
            {"denoised_dataset_address": str}.

        Raises:
            ValueError: If required settings are missing.
            requests.exceptions.RequestException: If the API request fails.
        """
        profile, project = self.validate_common_params(profile_name, project_name)

        data: Dict[str, Any] = {
            "profile_name": profile,
            "project_name": project,
            "dataset_address": dataset_address,
            "output_key": output_key,
            "strategy": strategy,
            "add_hit_labels": add_hit_labels,
            "hit_percentile": hit_percentile,
            "alpha": alpha,
            "drop_duplicates": drop_duplicates,
            "use_disynthon_pairs": use_disynthon_pairs,
            "aggregate_operation": aggregate_operation,
            "min_count_threshold": min_count_threshold,
        }

        if control_cols is not None:
            data["control_cols"] = control_cols
        if target_cols is not None:
            data["target_cols"] = target_cols
        if smiles_cols is not None:
            data["smiles_cols"] = smiles_cols

        response = self._post("/primitive/del/denoise", json=data)
        return self._validate_response(response)
