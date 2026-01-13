"""
Relative Binding Free Energy primitive module for DeepChem Server.

Contains the RelativeBindingFreeEnergyPrimitive class for relative binding free energy tasks.
"""

from typing import Any, Dict, Optional
from .base import Primitive


class RelativeBindingFreeEnergy(Primitive):
    """
    Primitive for relative binding free energy tasks.

    This class handles submitting relative binding free energy jobs to the DeepChem Server API.
    """

    def run(
        self,
        ligands_sdf_address: str,
        cleaned_protein_pdb_address: str,
        solvent: Optional[Dict] = None,
        overridden_rbfe_settings: Optional[Dict] = None,
        radial_network_central_ligand: Optional[str] = None,
        dry_run: bool = False,
        run_edges_in_parallel: bool = False,
        output_key: Optional[str] = "output_key",
        network_type: Optional[str] = "MINIMAL_SPANNING",
        scorer_type: Optional[str] = "LOMAP",
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Run the relative binding free energy primitive.

        Args:
            solvent: Dictionary of solvent settings
            ligands_sdf_address: Datastore address of ligands SDF file
            cleaned_protein_pdb_address: Datastore address of cleaned protein PDB file
            overridden_rbfe_settings: Dictionary of overridden RBFE settings
            radial_network_central_ligand: Central ligand for radial network
            dry_run: Whether to dry run the RBFE calculation
            run_edges_in_parallel: Whether to run edges in parallel
            output_key: Key for the output
            network_type: Type of network to use
            scorer_type: Type of scorer to use

        Returns:
            Response containing the relative binding free energy result

        Raises:
            ValueError: If required settings are missing
            requests.exceptions.RequestException: If API request fails
        """

        profile, project = self.validate_common_params(profile_name, project_name)

        data = {
            "profile_name": profile,
            "project_name": project,
            "ligands_sdf_address": ligands_sdf_address,
            "cleaned_protein_pdb_address": cleaned_protein_pdb_address,
            "solvent": solvent,
            "overridden_rbfe_settings": overridden_rbfe_settings,
            "radial_network_central_ligand": radial_network_central_ligand,
            "dry_run": dry_run,
            "run_edges_in_parallel": run_edges_in_parallel,
            "output_key": output_key,
            "network_type": network_type,
            "scorer_type": scorer_type,
        }

        response = self._post("/primitive/fep/calculate_rbfe", json=data, headers={"Content-Type": "application/json"})
        return self._validate_response(response)
