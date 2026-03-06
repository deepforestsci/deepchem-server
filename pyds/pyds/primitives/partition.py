"""
Partition primitive module for DeepChem Server.

Contains the Partition class for dataset partitioning tasks.
"""

from typing import Any, Dict, Optional

from .base import Primitive


class Partition(Primitive):
    """
    Primitive for dataset partitioning tasks.

    This class handles submitting partition jobs to the DeepChem Server API.
    """

    def run(
        self,
        dataset_address: str,
        n_partition: int = 4,
        shuffle: bool = False,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Run the partition primitive.

        Args:
            dataset_address: Datastore address of the dataset to partition
            n_partition: Number of partitions to create (default: 4)
            shuffle: Whether to shuffle before partitioning (only supported for DiskDataset)
            profile_name: Profile name (uses settings if not provided)
            project_name: Project name (uses settings if not provided)

        Returns:
            Response containing the partitioned dataset addresses

        Raises:
            ValueError: If required settings are missing
            requests.exceptions.RequestException: If API request fails
        """
        profile, project = self.validate_common_params(profile_name, project_name)

        data = {
            "profile_name": profile,
            "project_name": project,
            "dataset_address": dataset_address,
            "n_partition": n_partition,
            "shuffle": shuffle,
        }

        response = self._post("/primitive/partition", json=data)
        return self._validate_response(response)
