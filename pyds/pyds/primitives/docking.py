"""
Docking primitive module for DeepChem Server.

Contains the Docking class for VINA pose generation tasks.
"""

from typing import Optional

from ..models import Job
from .base import Primitive


class Docking(Primitive):
    """
    Primitive for molecular docking pose generation.

    This class submits docking jobs to the DeepChem Server API's
    /primitive/generate_pose endpoint.
    """

    def run(
        self,
        protein_address: str,
        ligand_address: str,
        output: str,
        exhaustiveness: int = 10,
        num_modes: int = 9,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Job:
        """
        Run the docking primitive.

        Args:
            protein_address: Datastore address of the protein PDB file
            ligand_address: Datastore address of the ligand file (PDB/SDF)
            output: Output name for docking results
            exhaustiveness: Vina exhaustiveness parameter (default: 10)
            num_modes: Number of binding modes to generate (default: 9)
            profile_name: Profile name (uses settings if not provided)
            project_name: Project name (uses settings if not provided)

        Returns:
            Job object representing the submitted docking job

        Raises:
            ValueError: If required settings are missing
            requests.exceptions.RequestException: If API request fails
        """
        profile, project = self.validate_common_params(profile_name, project_name)

        data = {
            "profile_name": profile,
            "project_name": project,
            "protein_address": protein_address,
            "ligand_address": ligand_address,
            "output": output,
            "exhaustiveness": exhaustiveness,
            "num_modes": num_modes,
        }

        response = self._post("/v1/primitive/generate_pose", json=data)
        response_data = self._validate_response(response)
        job = Job.from_dict(response_data, self)
        job.profile = profile
        job.project = project
        job.program_name = "generate_pose"
        return job
