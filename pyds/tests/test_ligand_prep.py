"""
Unit tests for LigandPrep primitive.
"""

import json

import pytest
import responses

from pyds.primitives import LigandPrep
from pyds.settings import Settings


class TestLigandPrep:
    """Unit tests for LigandPrep primitive."""

    def test_init(self, test_settings: Settings) -> None:
        """Test LigandPrep initialization."""
        client = LigandPrep(settings=test_settings)
        assert client.settings == test_settings
        assert client.base_url == "http://localhost:8000"

    @responses.activate
    def test_run_posts_expected_payload(self, test_settings: Settings) -> None:
        """Test POST body and endpoint for ligand-prep."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/ligand-prep",
            json={"ligand_sdf_address": "test_profile/test_project/out.sdf"},
            status=200,
        )

        client = LigandPrep(settings=test_settings)
        result = client.run(
            smiles="CCO",
            output="ethanol",
            ligand_name="ethanol",
            random_seed=7,
            optimize=False,
        )

        assert result["ligand_sdf_address"] == "test_profile/test_project/out.sdf"
        assert len(responses.calls) == 1
        body = json.loads(responses.calls[0].request.body)  # type: ignore
        assert body["profile_name"] == "test_profile"
        assert body["project_name"] == "test_project"
        assert body["smiles"] == "CCO"
        assert body["output"] == "ethanol"
        assert body["ligand_name"] == "ethanol"
        assert body["random_seed"] == 7
        assert body["optimize"] is False

    def test_run_missing_settings(self, test_settings_not_configured: Settings) -> None:
        """Test ligand_prep run with missing settings."""
        client = LigandPrep(settings=test_settings_not_configured)
        with pytest.raises(ValueError, match="Missing required settings"):
            client.run(smiles="CCO", output="out")
