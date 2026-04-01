"""
PDB cleaning primitive module for DeepChem Server.

Contains the PdbClean class for PDBFixer-based structure preparation.
"""

from typing import Any, Dict, List, Optional

from .base import Primitive


class PdbClean(Primitive):
    """
    Primitive for cleaning PDB structures via the DeepChem Server API.

    Submits jobs to the /primitive/pdb-clean endpoint.
    """

    def run(
        self,
        pdb_address: str,
        output: str,
        remove_chains: Optional[List[str]] = None,
        remove_heterogens: bool = True,
        remove_water: bool = True,
        add_hydrogens: bool = True,
        ph: float = 7.0,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Run the PDB cleaning primitive.

        Parameters
        ----------
        pdb_address : str
            Datastore address of the input PDB file.
        output : str
            Output key for the cleaned PDB in the datastore.
        remove_chains : list of str or None
            Chain IDs to remove before cleaning.
        remove_heterogens : bool
            Whether to strip heteroatom records.
        remove_water : bool
            Whether to remove water when stripping heterogens.
        add_hydrogens : bool
            Whether to add missing hydrogens.
        ph : float
            pH for protonation when adding hydrogens.
        profile_name : str or None
            Profile name; uses settings when omitted.
        project_name : str or None
            Project name; uses settings when omitted.

        Returns
        -------
        dict
            Response containing the cleaned_pdb_address entry.

        Raises
        ------
        ValueError
            If required settings are missing.
        """
        profile, project = self.validate_common_params(profile_name, project_name)

        data: Dict[str, Any] = {
            "profile_name": profile,
            "project_name": project,
            "pdb_address": pdb_address,
            "output": output,
            "remove_chains": remove_chains,
            "remove_heterogens": remove_heterogens,
            "remove_water": remove_water,
            "add_hydrogens": add_hydrogens,
            "ph": ph,
        }

        response = self._post("/primitive/pdb-clean", json=data)
        return self._validate_response(response)
