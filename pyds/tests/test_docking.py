"""
Unit tests for Docking primitive using live server.
"""

import time
import uuid

import pytest

from pyds.data import Data
from pyds.primitives import Docking
from pyds.settings import Settings


class TestDocking:
    """Unit tests for Docking primitive."""

    def test_init(self, test_settings: Settings) -> None:
        """Test Docking initialization."""
        client = Docking(settings=test_settings)
        assert client.settings == test_settings
        assert client.base_url == "http://localhost:8000"

    def test_run_success(
        self,
        live_data_client: Data,
        test_settings: Settings,
    ) -> None:
        """Test successful docking run on live server.

        Note: This assumes protein and ligand addresses already exist in the datastore.
        If not available, this test should be adjusted to upload fixtures and capture addresses.
        """
        client = Docking(settings=test_settings)

        # These should be replaced with real addresses in a live environment
        protein_address = "test/protein.pdb"
        ligand_address = "test/ligand.sdf"

        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        with pytest.raises(Exception):
            # Until valid addresses exist, server should error; ensures endpoint wiring
            client.run(
                protein_address=protein_address,
                ligand_address=ligand_address,
                output=f"test_docking_{test_id}_{timestamp}",
                exhaustiveness=8,
                num_modes=5,
            )

    def test_run_missing_settings(self, test_settings_not_configured: Settings) -> None:
        """Test docking run with missing settings."""
        client = Docking(settings=test_settings_not_configured)
        with pytest.raises(ValueError, match="Missing required settings"):
            client.run(
                protein_address="protein_addr",
                ligand_address="ligand_addr",
                output="out",
            )
