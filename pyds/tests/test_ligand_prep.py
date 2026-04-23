"""
Tests for LigandPrep primitive.
"""

import json
import os
import time
import uuid

import pytest
import responses
import tempfile
from pyds.data import Data
from pyds.primitives import LigandPrep
from pyds.settings import Settings


CORE_ASSETS = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    os.pardir,
    os.pardir,
    "deepchem_server",
    "core",
    "tests",
    "assets",
)

with open(os.path.join(CORE_ASSETS, "ligand_SMILES.json")) as _f:
    _LIGAND_DATA = json.load(_f)

SMILES_1 = "CC(=O)Oc1ccccc1C(=O)O"
SMILES_2 = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
SMILES_3 = "CCO"

PROTEOME_SCAN_SUBSET = {
    "Erlotinib": _LIGAND_DATA["ligand_smiles"]["Erlotinib"],
    "Selinexor": _LIGAND_DATA["ligand_smiles"]["Selinexor"],
    "Palbociclib": _LIGAND_DATA["ligand_smiles"]["Palbociclib"],
    "Olaparib": _LIGAND_DATA["ligand_smiles"]["Olaparib"],
}


def _unique_name(prefix: str) -> str:
    """Return a unique output name for each test invocation."""
    short_id = str(uuid.uuid4())[:8]
    return f"{prefix}_{short_id}_{int(time.time())}"


class TestLigandPrepUnit:
    """Unit tests for LigandPrep primitive (no server required)."""

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
            "http://localhost:8000/primitive/proteome-scan/ligand-prep",
            json={"ligand_sdf_address": "test_profile/test_project/out.sdf"},
            status=200,
        )

        client = LigandPrep(settings=test_settings)
        result = client.run(smiles="CCO", output="smiles_3", ligand_name="smiles_3", random_seed=7)

        assert result["ligand_sdf_address"] == "test_profile/test_project/out.sdf"
        assert len(responses.calls) == 1
        body = json.loads(responses.calls[0].request.body)  # type: ignore
        assert body["profile_name"] == "test_profile"
        assert body["project_name"] == "test_project"
        assert body["smiles"] == "CCO"
        assert body["output"] == "smiles_3"
        assert body["ligand_name"] == "smiles_3"
        assert body["random_seed"] == 7

    def test_run_missing_settings(self, test_settings_not_configured: Settings) -> None:
        """Test ligand_prep run with missing settings."""
        client = LigandPrep(settings=test_settings_not_configured)
        with pytest.raises(ValueError, match="Missing required settings"):
            client.run(smiles="CCO", output="out")


class TestLigandPrepLive:
    """Integration tests for LigandPrep using a live server."""

    def test_ligand_prep_basic(self, live_ligand_prep_client: LigandPrep) -> None:
        """Basic SDF generation from a valid SMILES."""
        result = live_ligand_prep_client.run(
            smiles=SMILES_1,
            output=_unique_name("smiles_1"),
        )
        assert "ligand_sdf_address" in result

    def test_ligand_prep_output_extension(self, live_ligand_prep_client: LigandPrep) -> None:
        """Output key gets .sdf extension when not supplied."""
        result = live_ligand_prep_client.run(
            smiles=SMILES_3,
            output=_unique_name("smiles_3"),
        )
        assert result["ligand_sdf_address"].endswith(".sdf")

        result2 = live_ligand_prep_client.run(
            smiles=SMILES_3,
            output=_unique_name("smiles_32") + ".sdf",
        )
        assert result2["ligand_sdf_address"].endswith(".sdf")

    def test_ligand_prep_valid_sdf_content(
        self,
        live_ligand_prep_client: LigandPrep,
        live_data_client: Data,
    ) -> None:
        """Retrieved SDF file contains expected molecule keywords."""
        result = live_ligand_prep_client.run(
            smiles=SMILES_1,
            output=_unique_name("smiles_1_content"),
        )

        addr = result["ligand_sdf_address"]
        filename = addr.split("/")[-1] if "/" in addr else addr

        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            live_data_client.download_data(address=filename, destination_path=tmp_path)
            with open(tmp_path, "r") as f:
                sdf_text = f.read()
            assert "M  END" in sdf_text
            assert "$$$$" in sdf_text
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def test_ligand_prep_with_ligand_name(
        self,
        live_ligand_prep_client: LigandPrep,
        live_data_client: Data,
    ) -> None:
        """Ligand name is embedded in the SDF block."""
        result = live_ligand_prep_client.run(
            smiles=SMILES_1,
            output=_unique_name("smiles_1_named"),
            ligand_name="smiles_1",
        )

        addr = result["ligand_sdf_address"]
        filename = addr.split("/")[-1] if "/" in addr else addr

        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            live_data_client.download_data(address=filename, destination_path=tmp_path)
            with open(tmp_path, "r") as f:
                sdf_text = f.read()
            assert "smiles_1" in sdf_text
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def test_ligand_prep_without_seed(self, live_ligand_prep_client: LigandPrep) -> None:
        """Primitive succeeds without specifying a random seed (uses default)."""
        result = live_ligand_prep_client.run(
            smiles=SMILES_1,
            output=_unique_name("smiles_1_noseed"),
        )
        assert "ligand_sdf_address" in result

    def test_ligand_prep_reproducible_with_seed(
        self,
        live_ligand_prep_client: LigandPrep,
    ) -> None:
        """Same SMILES and seed produces valid results."""
        result1 = live_ligand_prep_client.run(
            smiles=SMILES_1,
            output=_unique_name("smiles_1_seed1"),
            random_seed=42,
        )
        result2 = live_ligand_prep_client.run(
            smiles=SMILES_1,
            output=_unique_name("smiles_1_seed2"),
            random_seed=42,
        )
        assert "ligand_sdf_address" in result1
        assert "ligand_sdf_address" in result2

    def test_ligand_prep_invalid_smiles(self, live_ligand_prep_client: LigandPrep) -> None:
        """Invalid SMILES causes the server to return an error."""
        with pytest.raises(Exception):
            live_ligand_prep_client.run(
                smiles="not_a_smiles!!!",
                output=_unique_name("bad_mol"),
            )

    def test_ligand_prep_empty_smiles(self, live_ligand_prep_client: LigandPrep) -> None:
        """Empty SMILES causes the server to return an error."""
        with pytest.raises(Exception):
            live_ligand_prep_client.run(smiles="", output=_unique_name("empty"))

    def test_ligand_prep_different_molecules(
        self,
        live_ligand_prep_client: LigandPrep,
    ) -> None:
        """Multiple different molecules all produce valid results."""
        for name, smiles in [
            ("smiles1", SMILES_1),
            ("smiles2", SMILES_2),
            ("smiles3", SMILES_3),
        ]:
            result = live_ligand_prep_client.run(
                smiles=smiles,
                output=_unique_name(f"mol_{name}"),
            )
            assert "ligand_sdf_address" in result
            assert result["ligand_sdf_address"].endswith(".sdf")

    def test_ligand_prep_proteome_scan_drugs(
        self,
        live_ligand_prep_client: LigandPrep,
    ) -> None:
        """Subset of 4 ProteomeScan clinical drug molecules all succeed."""
        for name, smiles in PROTEOME_SCAN_SUBSET.items():
            result = live_ligand_prep_client.run(
                smiles=smiles,
                output=_unique_name(f"ps_{name.lower()}"),
                ligand_name=name,
            )
            assert "ligand_sdf_address" in result
            assert result["ligand_sdf_address"].endswith(".sdf")

    def test_ligand_prep_drug_names_in_sdf(
        self,
        live_ligand_prep_client: LigandPrep,
        live_data_client: Data,
    ) -> None:
        """Ligand names from the ProteomeScan dataset appear in the SDF block."""

        name = "Erlotinib"
        smiles = PROTEOME_SCAN_SUBSET[name]
        result = live_ligand_prep_client.run(
            smiles=smiles,
            output=_unique_name(f"named_{name.lower()}"),
            ligand_name=name,
        )

        addr = result["ligand_sdf_address"]
        filename = addr.split("/")[-1] if "/" in addr else addr

        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            live_data_client.download_data(address=filename, destination_path=tmp_path)
            with open(tmp_path, "r") as f:
                sdf_text = f.read()
            assert name in sdf_text, f"Ligand name '{name}' not embedded in SDF"
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
