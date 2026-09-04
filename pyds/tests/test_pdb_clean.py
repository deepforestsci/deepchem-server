"""
Tests for PdbClean primitive.
"""

import os
import tempfile
import time
import uuid

import pytest
import responses
from pyds.data import Data
from pyds.primitives import PdbClean
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
RAW_PDB = os.path.join(CORE_ASSETS, "cleaned_3cyx.pdb")
TWO_CHAIN_PDB = os.path.join(CORE_ASSETS, "test_protein_2chain.pdb")


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
        """Test POST body and endpoint for pdb-clean."""
        import json

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


def _upload_pdb(data_client: Data, pdb_path: str, filename: str) -> str:
    """Upload a PDB file using the Data client and return the deepchem address."""
    result = data_client.upload_data(file_path=pdb_path, filename=filename)
    return result["dataset_address"]


class TestPdbCleanLive:
    """Integration tests for PdbClean using a live server."""

    def test_pdb_clean_basic(
        self,
        live_pdb_clean_client: PdbClean,
        live_data_client: Data,
    ) -> None:
        """Basic PDB cleaning returns a cleaned_pdb_address."""
        pdb_addr = _upload_pdb(live_data_client, RAW_PDB, "raw_protein.pdb")
        result = live_pdb_clean_client.run(
            pdb_address=pdb_addr,
            output=_unique_name("cleaned_basic"),
        )
        addr = result["cleaned_pdb_address"]
        assert "cleaned_basic" in addr

    def test_pdb_clean_output_extension(
        self,
        live_pdb_clean_client: PdbClean,
        live_data_client: Data,
    ) -> None:
        """Output key gets .pdb extension when not supplied."""
        pdb_addr = _upload_pdb(live_data_client, RAW_PDB, "raw_ext_test.pdb")

        result_no_ext = live_pdb_clean_client.run(
            pdb_address=pdb_addr,
            output=_unique_name("cleaned_no_ext"),
        )
        assert result_no_ext["cleaned_pdb_address"].endswith(".pdb")

        result_with_ext = live_pdb_clean_client.run(
            pdb_address=pdb_addr,
            output=_unique_name("cleaned_with_ext") + ".pdb",
        )
        assert result_with_ext["cleaned_pdb_address"].endswith(".pdb")

    def test_pdb_clean_produces_valid_pdb(
        self,
        live_pdb_clean_client: PdbClean,
        live_data_client: Data,
    ) -> None:
        """Cleaned PDB file contains ATOM records when downloaded."""
        pdb_addr = _upload_pdb(live_data_client, RAW_PDB, "raw_valid_test.pdb")
        result = live_pdb_clean_client.run(
            pdb_address=pdb_addr,
            output=_unique_name("cleaned_valid"),
        )

        cleaned_addr = result["cleaned_pdb_address"]
        filename = cleaned_addr.split("/")[-1] if "/" in cleaned_addr else cleaned_addr
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            live_data_client.download_data(address=filename, destination_path=tmp_path)
            with open(tmp_path, "r") as f:
                content = f.read()
            assert "ATOM" in content or "HETATM" in content
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def test_pdb_clean_no_heterogens(
        self,
        live_pdb_clean_client: PdbClean,
        live_data_client: Data,
    ) -> None:
        """Cleaned PDB with remove_heterogens=True and remove_water=True has no HOH."""
        pdb_addr = _upload_pdb(live_data_client, RAW_PDB, "raw_no_hets.pdb")
        result = live_pdb_clean_client.run(
            pdb_address=pdb_addr,
            output=_unique_name("cleaned_no_hets"),
            remove_heterogens=True,
            remove_water=True,
        )

        filename = result["cleaned_pdb_address"].split("/")[-1]
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            live_data_client.download_data(address=filename, destination_path=tmp_path)
            with open(tmp_path, "r") as f:
                content = f.read()
            water_lines = [line for line in content.splitlines() if line.startswith("HETATM") and "HOH" in line]
            assert len(water_lines) == 0
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def test_pdb_clean_custom_ph(
        self,
        live_pdb_clean_client: PdbClean,
        live_data_client: Data,
    ) -> None:
        """Cleaning with a non-default pH completes without error."""
        pdb_addr = _upload_pdb(live_data_client, RAW_PDB, "raw_ph_test.pdb")
        result = live_pdb_clean_client.run(
            pdb_address=pdb_addr,
            output=_unique_name("cleaned_ph6"),
            ph=6.0,
        )
        assert "cleaned_pdb_address" in result

    def test_pdb_clean_without_hydrogens(
        self,
        live_pdb_clean_client: PdbClean,
        live_data_client: Data,
    ) -> None:
        """add_hydrogens=False still produces a valid result."""
        pdb_addr = _upload_pdb(live_data_client, RAW_PDB, "raw_noh_test.pdb")
        result = live_pdb_clean_client.run(
            pdb_address=pdb_addr,
            output=_unique_name("cleaned_noh"),
            add_hydrogens=False,
        )
        assert "cleaned_pdb_address" in result

    def test_pdb_clean_empty_address_raises(
        self,
        live_pdb_clean_client: PdbClean,
    ) -> None:
        """Empty pdb_address causes the server to return an error."""
        with pytest.raises(Exception):
            live_pdb_clean_client.run(pdb_address="", output="should_fail")

    def test_pdb_clean_two_chain_protein(
        self,
        live_pdb_clean_client: PdbClean,
        live_data_client: Data,
    ) -> None:
        """Cleaning the two-chain fixture succeeds."""
        pdb_addr = _upload_pdb(live_data_client, TWO_CHAIN_PDB, "two_chain.pdb")
        result = live_pdb_clean_client.run(
            pdb_address=pdb_addr,
            output=_unique_name("cleaned_two_chain"),
        )
        assert "cleaned_pdb_address" in result

    def test_pdb_clean_two_chain_has_atom_records(
        self,
        live_pdb_clean_client: PdbClean,
        live_data_client: Data,
    ) -> None:
        """Cleaned two-chain structure contains ATOM records."""
        pdb_addr = _upload_pdb(live_data_client, TWO_CHAIN_PDB, "two_chain_atoms.pdb")
        result = live_pdb_clean_client.run(
            pdb_address=pdb_addr,
            output=_unique_name("cleaned_two_chain_atoms"),
        )

        filename = result["cleaned_pdb_address"].split("/")[-1]
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            live_data_client.download_data(address=filename, destination_path=tmp_path)
            with open(tmp_path, "r") as f:
                content = f.read()
            atom_lines = [line for line in content.splitlines() if line.startswith("ATOM")]
            assert len(atom_lines) > 0
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def test_pdb_clean_water_retained(
        self,
        live_pdb_clean_client: PdbClean,
        live_data_client: Data,
    ) -> None:
        """remove_water=False preserves HOH records."""
        pdb_addr = _upload_pdb(live_data_client, TWO_CHAIN_PDB, "two_chain_keepwater.pdb")
        result = live_pdb_clean_client.run(
            pdb_address=pdb_addr,
            output=_unique_name("cleaned_keepwater"),
            remove_heterogens=True,
            remove_water=False,
        )

        filename = result["cleaned_pdb_address"].split("/")[-1]
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            live_data_client.download_data(address=filename, destination_path=tmp_path)
            with open(tmp_path, "r") as f:
                content = f.read()
            hoh_lines = [line for line in content.splitlines() if "HOH" in line]
            assert len(hoh_lines) > 0, "Expected HOH water records to be retained"
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def test_pdb_clean_heterogens_kept_when_disabled(
        self,
        live_pdb_clean_client: PdbClean,
        live_data_client: Data,
    ) -> None:
        """remove_heterogens=False leaves all HETATM records intact."""
        pdb_addr = _upload_pdb(live_data_client, TWO_CHAIN_PDB, "two_chain_keephets.pdb")
        result = live_pdb_clean_client.run(
            pdb_address=pdb_addr,
            output=_unique_name("cleaned_keephets"),
            remove_heterogens=False,
            remove_water=False,
        )

        filename = result["cleaned_pdb_address"].split("/")[-1]
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            live_data_client.download_data(address=filename, destination_path=tmp_path)
            with open(tmp_path, "r") as f:
                content = f.read()
            hetatm_lines = [line for line in content.splitlines() if line.startswith("HETATM")]
            assert len(hetatm_lines) > 0, "Expected HETATM records to be kept"
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
