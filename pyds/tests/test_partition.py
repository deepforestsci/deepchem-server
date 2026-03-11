"""
Unit tests for Partition primitive using live server.
"""

import time
import uuid

import pytest

from pyds.data import Data
from pyds.primitives import Partition
from pyds.settings import Settings


class TestPartition:
    """Unit tests for Partition primitive."""

    def test_init(self, test_settings: Settings) -> None:
        """Test Partition initialization."""
        client = Partition(settings=test_settings)

        assert client.settings == test_settings
        assert client.base_url == "http://localhost:8000"

    def test_run_success(
        self,
        live_partition_client: Partition,
        live_data_client: Data,
        small_classification_csv: str,
    ) -> None:
        """Test successful partition run on live server."""
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        upload_result = live_data_client.upload_data(
            file_path=small_classification_csv,
            filename=f"test_partition_{test_id}_{timestamp}.csv",
            description="Test data for partitioning",
        )
        dataset_address = upload_result["dataset_address"]

        result = live_partition_client.run(
            dataset_address=dataset_address,
            n_partition=2,
        )

        assert "partitioned_dataset_addresses" in result
        assert len(result["partitioned_dataset_addresses"]) == 2

    def test_run_with_defaults(
        self,
        live_partition_client: Partition,
        live_data_client: Data,
        small_classification_csv: str,
    ) -> None:
        """Test partition run with default parameters on live server."""
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        upload_result = live_data_client.upload_data(
            file_path=small_classification_csv,
            filename=f"test_partition_defaults_{test_id}_{timestamp}.csv",
            description="Test data for partitioning with defaults",
        )
        dataset_address = upload_result["dataset_address"]

        result = live_partition_client.run(dataset_address=dataset_address,)

        assert "partitioned_dataset_addresses" in result

    def test_run_with_profile_project_override(
        self,
        live_partition_client: Partition,
        live_data_client: Data,
        small_classification_csv: str,
    ) -> None:
        """Test partition run with profile and project override on live server."""
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        upload_result = live_data_client.upload_data(
            file_path=small_classification_csv,
            filename=f"test_partition_override_{test_id}_{timestamp}.csv",
            description="Test data for partitioning with override",
        )
        dataset_address = upload_result["dataset_address"]

        result = live_partition_client.run(
            dataset_address=dataset_address,
            n_partition=2,
            profile_name="test_profile",
            project_name="test_project",
        )

        assert "partitioned_dataset_addresses" in result

    def test_run_missing_settings(self, test_settings_not_configured: Settings) -> None:
        """Test partition run with missing settings."""
        client = Partition(settings=test_settings_not_configured)

        with pytest.raises(ValueError, match="Missing required settings"):
            client.run(
                dataset_address="test/dataset",
                n_partition=4,
            )

    def test_run_api_error(self, live_partition_client: Partition) -> None:
        """Test partition run with API error on live server."""
        with pytest.raises(Exception):
            live_partition_client.run(
                dataset_address="non_existent_dataset",
                n_partition=4,
            )
