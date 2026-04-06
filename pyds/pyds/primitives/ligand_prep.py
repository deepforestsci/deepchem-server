"""
Ligand preparation primitive module for DeepChem Server.

Contains the LigandPrep class for SMILES to 3D SDF conversion jobs.
"""

from typing import Any, Dict, Optional

from .base import Primitive


class LigandPrep(Primitive):
    """
    Primitive for preparing ligands (SMILES to 3D SDF) via the DeepChem Server API.

    Submits jobs to the /primitive/ligand-prep endpoint.
    """

    def run(
        self,
        smiles: str,
        output: str,
        ligand_name: str = "",
        random_seed: int = 42,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Run the ligand preparation primitive.

        Parameters
        ----------
        smiles : str
            Input SMILES string.
        output : str
            Output key for the SDF in the datastore.
        ligand_name : str
            Optional molecule name stored in the SDF.
        random_seed : int
            Seed for 3D embedding.
        optimize : bool
            Whether to run MMFF94 optimization.
        profile_name : str or None
            Profile name; uses settings when omitted.
        project_name : str or None
            Project name; uses settings when omitted.

        Returns
        -------
        dict
            Response containing the ligand_sdf_address entry.

        Raises
        ------
        ValueError
            If required settings are missing.
        """
        profile, project = self.validate_common_params(profile_name, project_name)

        data: Dict[str, Any] = {
            "profile_name": profile,
            "project_name": project,
            "smiles": smiles,
            "output": output,
            "ligand_name": ligand_name,
            "random_seed": random_seed
        }

        response = self._post("/primitive/ligand-prep", json=data)
        return self._validate_response(response)
