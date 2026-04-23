"""
Integration tests for the ProteomeScan client.

These tests verify the high-level client against a running deepchem
server. They exercise:

1. Client initialization
2. Individual primitive helpers (``get_gene_pdb``, ``dock_pair``)
3. The full ``run_pipeline`` orchestrator

All tests require a running deepchem-server (via ``live_settings``).
"""

import json
import os
import tempfile
from io import StringIO

import pandas as pd
import pytest

from pyds.data import Data
from pyds.primitives.proteome_scan import ProteomeScan


@pytest.fixture
def test_ligand_sdf():
    """Create a test ligand SDF file."""
    sdf_content = (
        "test_ligand\n"
        "  test\n\n"
        "  6  6  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.4000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    2.1000    1.2124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.4000    2.4249    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.0000    2.4249    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.7000    1.2124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  2  3  1  0  0  0  0\n"
        "  3  4  1  0  0  0  0\n"
        "  4  5  1  0  0  0  0\n"
        "  5  6  1  0  0  0  0\n"
        "  6  1  1  0  0  0  0\n"
        "M  END\n"
        "$$$$\n")
    with tempfile.NamedTemporaryFile(
            mode="w", suffix=".sdf", delete=False) as f:
        f.write(sdf_content)
        temp_file = f.name

    yield temp_file

    if os.path.exists(temp_file):
        os.unlink(temp_file)


@pytest.fixture
def uploaded_ligand(live_data_client, test_ligand_sdf):
    """Upload the test ligand to the live datastore."""
    result = live_data_client.upload_data(test_ligand_sdf)
    return result["dataset_address"]


class TestProteomeScanClient:
    """Integration tests for the ProteomeScan client."""

    def test_client_initialization(self, live_settings):
        client = ProteomeScan(settings=live_settings)
        assert client.settings is live_settings
        assert client.pdb_clean is not None
        assert client.run_docking is not None
        assert client.parse_results is not None
        assert client.run_multi_pose_analysis is not None

    def test_get_gene_pdb(self, live_proteome_scan_client):
        gene_name = "MB"
        entry_id = "P02144"
        scan_id = "test_client_get_gene_pdb"

        result_addr = live_proteome_scan_client.get_gene_pdb(
            gene_name=gene_name,
            entry_id=entry_id,
            scan_id=scan_id,
        )

        assert result_addr.startswith("deepchem://")

        data_client = Data(live_proteome_scan_client.settings)
        result_data = data_client.get_data(result_addr)
        if isinstance(result_data, str):
            result_data = json.loads(result_data)

        assert result_data["gene_name"] == gene_name
        assert result_data["scan_id"] == scan_id
        assert "pdbs_csv_address" in result_data

    def test_dock_pair(self, live_proteome_scan_client, uploaded_ligand):
        gene_name = "MB"
        entry_id = "P02144"
        scan_id = "test_client_dock_pair"
        ligand_name = "TestLigand"

        live_proteome_scan_client.get_gene_pdb(
            gene_name=gene_name,
            entry_id=entry_id,
            scan_id=scan_id,
        )

        result_addr = live_proteome_scan_client.dock_pair(
            gene_name=gene_name,
            ligand_name=ligand_name,
            ligand_address=uploaded_ligand,
            scan_id=scan_id,
            exhaustiveness=1,
            num_modes=1,
        )

        assert result_addr.startswith("deepchem://")
        data_client = Data(live_proteome_scan_client.settings)
        result_data = data_client.get_data(result_addr)
        if isinstance(result_data, str):
            result_data = json.loads(result_data)
        assert result_data["gene_name"] == gene_name
        assert result_data["ligand_name"] == ligand_name
        assert "top_score_csv_address" in result_data

    def test_run_pipeline_basic(self, live_proteome_scan_client,
                                uploaded_ligand):
        scan_id = "test_client_pipeline"
        ligands = ["TestLigand"]
        gene_names = ["MB"]
        gene_entry_map = {"MB": "P02144"}
        ligand_address_map = {"TestLigand": uploaded_ligand}

        result = live_proteome_scan_client.run_pipeline(
            scan_id=scan_id,
            ligands=ligands,
            gene_names=gene_names,
            gene_entry_map=gene_entry_map,
            ligand_address_map=ligand_address_map,
            exhaustiveness=1,
            num_modes=1,
        )

        assert result["scan_id"] == scan_id
        assert "gene_pdbs" in result
        assert "docking_results" in result
        assert "parsed_results" in result

        assert "MB" in result["gene_pdbs"]
        assert "MB" in result["docking_results"]

        data_client = Data(live_proteome_scan_client.settings)
        parsed_data = data_client.get_data(result["parsed_results"])
        if isinstance(parsed_data, str):
            parsed_data = json.loads(parsed_data)

        assert parsed_data["scan_id"] == scan_id
        assert "aggregated" in parsed_data
        agg = parsed_data["aggregated"]["TestLigand"]
        assert agg["address"].startswith("deepchem://")

        df = data_client.get_data(agg["address"])
        if isinstance(df, str):
            df = pd.read_csv(StringIO(df))
        assert len(df) > 0
        assert "gene_name" in df.columns
        assert "top_score" in df.columns

    def test_run_pipeline_invalid_gene(self, live_proteome_scan_client,
                                       uploaded_ligand):
        scan_id = "test_client_pipeline_invalid"
        ligands = ["TestLigand"]
        gene_names = ["INVALID_GENE"]
        gene_entry_map = {"INVALID_GENE": "INVALID123"}
        ligand_address_map = {"TestLigand": uploaded_ligand}

        with pytest.raises(Exception):
            live_proteome_scan_client.run_pipeline(
                scan_id=scan_id,
                ligands=ligands,
                gene_names=gene_names,
                gene_entry_map=gene_entry_map,
                ligand_address_map=ligand_address_map,
            )
