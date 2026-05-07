"""
ProteomeScan primitive wrappers for DeepChem Server.
"""

from typing import Any, Dict, Optional

from .base import Primitive


class PdbClean(Primitive):
    """Downloads and cleans the best experimental PDB structures for a gene.

    Submits a job to /primitive/pdb_clean. The server fetches all experimental
    structures for the given UniProt entry, selects the optimal non-overlapping
    coverage set, cleans each with PDBFixer/OpenMM, and returns a JSON summary
    with datastore addresses for all cleaned files.
    """

    def run(
        self,
        gene_name: str,
        entry_id: str,
        scan_id: str,
        output: str,
        min_res_val: float = 2.5,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Run the pdb_clean primitive.

        Parameters
        ----------
        gene_name : str
            Human gene name (e.g. 'MAP2K1').
        entry_id : str
            UniProt accession (e.g. 'P02144').
        scan_id : str
            Identifier for the current scan run.
        output : str
            Prefix used when naming files in the datastore.
        min_res_val : float, optional
            Resolution threshold (Angstrom) for the high-res candidate set.
            Default is 2.5.
        profile_name : str, optional
            Profile name (uses settings if not provided).
        project_name : str, optional
            Project name (uses settings if not provided).

        Returns
        -------
        Dict[str, Any]
            {"pdb_clean_results_address": str}

        Raises
        ------
        ValueError
            If required settings are missing.
        requests.exceptions.RequestException
            If the API request fails.

        Examples
        --------
        >>> client = PdbClean(settings)
        >>> result = client.run("MB", "P02144", "scan_001", "mb_output")
        >>> result["pdb_clean_results_address"]
        'deepchem://profile/project/mb_output/MB_pdb_clean_summary.json'
        """
        profile, project = self.validate_common_params(profile_name, project_name)

        data = {
            "profile_name": profile,
            "project_name": project,
            "gene_name": gene_name,
            "entry_id": entry_id,
            "scan_id": scan_id,
            "output": output,
            "min_res_val": min_res_val,
        }

        response = self._post("/primitive/pdb_clean", json=data)
        return self._validate_response(response)
