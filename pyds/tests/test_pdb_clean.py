"""
Unit tests for PdbClean primitive.
"""

import json

import pytest
import responses

from pyds.primitives import PdbClean
from pyds.settings import Settings


class TestPdbClean:
    """Unit tests for PdbClean primitive."""

    def test_init(self, test_settings: Settings) -> None:
        """Test PdbClean initialization."""
        client = PdbClean(settings=test_settings)
        assert client.settings == test_settings
        assert client.base_url == "http://localhost:8000"

    @responses.activate
    def test_run_posts_expected_payload(self, test_settings: Settings) -> None:
        """Test POST body and endpoint for pdb-clean."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/pdb-clean",
            json={"cleaned_pdb_address": "test_profile/test_project/cleaned.pdb"},
            status=200,
        )

        client = PdbClean(settings=test_settings)
        result = client.run(
            pdb_address="deepchem://test_profile/test_project/raw.pdb",
            output="cleaned",
            remove_chains=["B"],
            remove_heterogens=True,
            remove_water=False,
            add_hydrogens=True,
            ph=6.5,
        )

        assert result["cleaned_pdb_address"] == "test_profile/test_project/cleaned.pdb"
        assert len(responses.calls) == 1
        body = json.loads(responses.calls[0].request.body)  # type: ignore
        assert body["profile_name"] == "test_profile"
        assert body["project_name"] == "test_project"
        assert body["pdb_address"] == "deepchem://test_profile/test_project/raw.pdb"
        assert body["output"] == "cleaned"
        assert body["remove_chains"] == ["B"]
        assert body["remove_heterogens"] is True
        assert body["remove_water"] is False
        assert body["add_hydrogens"] is True
        assert body["ph"] == 6.5

    def test_run_missing_settings(self, test_settings_not_configured: Settings) -> None:
        """Test pdb_clean run with missing settings."""
        client = PdbClean(settings=test_settings_not_configured)
        with pytest.raises(ValueError, match="Missing required settings"):
            client.run(pdb_address="addr", output="out")
