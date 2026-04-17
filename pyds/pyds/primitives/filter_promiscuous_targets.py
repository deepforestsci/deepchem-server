"""
FilterPromiscuousTargets primitive module for DeepChem Server.
"""

from typing import Any, Dict, List, Optional

from .base import Primitive


class FilterPromiscuousTargets(Primitive):
    """
    Primitive for filtering promiscuous targets from docking scan results.

    This class submits jobs to the DeepChem Server API's
    /primitive/filter-promiscuous-targets endpoint.
    """

    def run(
        self,
        scan_result_addresses: List[str],
        thresholds: List[List[int]],
        output: str,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Run the filter_promiscuous_targets primitive.

        Parameters
        ----------
        scan_result_addresses: List[str]
            DeepChem addresses of per-ligand scan result
                CSV files. Each CSV must contain a gene_name column and rows
                should be sorted from best (top) to worst docking score.
        thresholds: List[List[int]]
            List of [m, n] pairs where m is the percentile cutoff
                (0-100) and n is the minimum occurrence count across ligands.
                For example, [[15, 1]] marks a gene as promiscuous if it
                appears in the top 15% of at least 1 ligand's results.
        output: str
            Output name prefix used for all uploaded result files.
        profile_name: Optional[str]
            profile_name: Profile name (uses settings if not provided).
        project_name: Optional[str]
            Project name (uses settings if not provided).

        Returns
        -------
            Response dict with key filter_results_address pointing to a
            JSON file that contains:

            - promiscuous_targets: dict mapping threshold string to list of
              promiscuous gene names
            - promiscuous_targets_address: address of standalone JSON
            - filtered_results: dict mapping threshold string to a dict of
              {original_filename: filtered_csv_address}

        Raises
        ------
            ValueError: If required settings are missing.
            requests.exceptions.RequestException: If API request fails.

        Examples
        --------
        >>> result = FilterPromiscuousTargets().run(
        ...     scan_result_addresses=["deepchem://test_profile/test_project/scan_a.csv", "deepchem://test_profile/test_project/scan_b.csv"],
        ...     thresholds=[[15, 1]],
        ...     output="test_output",
        ... )
        >>> print(result["filter_results_address"])
        deepchem://test_profile/test_project/test_output_results.json
        """
        profile, project = self.validate_common_params(profile_name, project_name)

        data = {
            "profile_name": profile,
            "project_name": project,
            "scan_result_addresses": scan_result_addresses,
            "thresholds": thresholds,
            "output": output,
        }

        response = self._post("/primitive/filter-promiscuous-targets", json=data)
        return self._validate_response(response)
