"""
Integration tests for the PdbClean primitive using a live server.
"""

import time
import uuid

import pytest

from pyds.primitives.proteome_scan import PdbClean
from pyds.settings import Settings
import requests


class TestPdbClean:
    """Integration tests for PdbClean primitive."""

    def test_init(self, test_settings: Settings) -> None:
        """Test PdbClean initialization."""
        client = PdbClean(settings=test_settings)
        assert client.settings == test_settings
        assert client.base_url == "http://localhost:8000"

    def test_pdb_clean_mb_p02144(
        self,
        live_pdb_clean_client: PdbClean,
    ) -> None:
        """Run pdb_clean against a real UniProt entry on a live server."""
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        try:
            result = live_pdb_clean_client.run(
                gene_name="MB",
                entry_id="P02144",
                scan_id=f"pyds_pdb_clean_{test_id}_{timestamp}",
                output=f"pyds_pdb_clean_out_{test_id}_{timestamp}",
                min_res_val=2.5,
            )
        except requests.exceptions.HTTPError as e:
            msg = str(e)
            if "PROTEOMESCAN_CACHE_ROOT environment variable is not set" in msg:
                pytest.skip("Server is missing PROTEOMESCAN_CACHE_ROOT. "
                            "Set it on the deepchem-server process to enable pdb_clean integration tests.")
            raise

        assert "pdb_clean_results_address" in result
        assert isinstance(result["pdb_clean_results_address"], str)
        assert len(result["pdb_clean_results_address"]) > 0

    def test_missing_settings(self, test_settings_not_configured: Settings) -> None:
        """Test pdb_clean run with missing settings."""
        client = PdbClean(settings=test_settings_not_configured)

        with pytest.raises(ValueError, match="Missing required settings"):
            client.run(
                gene_name="MB",
                entry_id="P02144",
                scan_id="scan",
                output="output",
            )
