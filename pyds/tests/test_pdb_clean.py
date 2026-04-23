"""
Tests for PdbClean primitive.
"""

import os
import time
import uuid

import pytest
import responses
from pyds.primitives import PdbClean
from pyds.settings import Settings


def _unique_name(prefix: str) -> str:
    """Return a unique output name for each test invocation."""
    short_id = str(uuid.uuid4())[:8]
    return f"{prefix}_{short_id}_{int(time.time())}"


class TestPdbCleanUnit:
    """Unit tests for PdbClean primitive (no server required)."""

    def test_init(self, test_settings: Settings) -> None:
        """Test PdbClean initialization."""
        client = PdbClean(settings=test_settings)
        assert client.settings == test_settings
        assert client.base_url == "http://localhost:8000"

    @responses.activate
    def test_run_posts_expected_payload(self, test_settings: Settings) -> None:
        """Test POST body and endpoint for proteome-scan/pdb-clean."""
        import json

        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/proteome-scan/pdb-clean",
            json={"results_address": "test_profile/test_project/pdb_clean/results.json"},
            status=200,
        )

        client = PdbClean(settings=test_settings)
        result = client.run(
            gene_name="GBA3",
            entry_id="Q9H3H0",
            scan_id="scan_test_001",
            output="gba3_run",
            min_res_val=2.0,
        )

        assert result["results_address"] == "test_profile/test_project/pdb_clean/results.json"
        assert len(responses.calls) == 1
        body = json.loads(responses.calls[0].request.body)  # type: ignore
        assert body["profile_name"] == "test_profile"
        assert body["project_name"] == "test_project"
        assert body["gene_name"] == "GBA3"
        assert body["entry_id"] == "Q9H3H0"
        assert body["scan_id"] == "scan_test_001"
        assert body["output"] == "gba3_run"
        assert body["min_res_val"] == 2.0

    def test_run_missing_settings(self, test_settings_not_configured: Settings) -> None:
        """Test pdb_clean run with missing settings."""
        client = PdbClean(settings=test_settings_not_configured)
        with pytest.raises(ValueError, match="Missing required settings"):
            client.run(gene_name="GBA3", entry_id="Q9H3H0", scan_id="scan_x", output="out")


class TestPdbCleanLive:
    """Integration tests for PdbClean using a live server."""

    def test_pdb_clean_basic(
        self,
        live_pdb_clean_client: PdbClean,
    ) -> None:
        """Basic PDB cleaning returns a results_address."""
        result = live_pdb_clean_client.run(
            # Q9H3H0 can have no PDBe/PDB mappings and may 404 on PDBe.
            # P69905 (human hemoglobin subunit alpha) has stable PDB coverage.
            gene_name="HBA1",
            entry_id="P69905",
            scan_id=_unique_name("scan"),
            output=_unique_name("hba1"),
        )
        assert "results_address" in result
