"""Integration tests for the pdb_clean primitive.

These tests do not mock external dependencies. They perform real network calls
to EMBL-EBI Proteins API and PDBe endpoints and run the real PDBFixer/OpenMM
cleaning step.
"""

from __future__ import annotations

import json
from pathlib import Path
import shutil

import pandas as pd
import pytest

from deepchem_server.core import config
from deepchem_server.core.primitives.proteome_scan.pdb_clean import (
    _Range,
    _get_optimal_pdbs_df,
    _select_ranges,
    pdb_clean,
)


ASSETS_DIR = Path(__file__).parent.parent / "assets" / "proteome_scan"
MINIMAL_PDB = ASSETS_DIR / "minimal_protein.pdb"


class TestSelectRanges:

    def test_non_overlapping_all_selected(self):
        ranges = [
            _Range("A", 1, 50, 2.0),
            _Range("B", 60, 120, 2.5),
            _Range("C", 130, 200, 1.8),
        ]
        result = _select_ranges(ranges)
        assert len(result) == 3

    def test_overlapping_better_resolution_replaces(self):
        ranges = [
            _Range("A", 1, 100, 3.0),
            _Range("B", 50, 150, 1.5),
        ]
        result = _select_ranges(ranges)
        # B has better resolution and more coverage so it should be selected
        ids = [r.id for r in result]
        assert "B" in ids

    def test_empty_input(self):
        assert _select_ranges([]) == []

    def test_single_range(self):
        r = _Range("A", 1, 100, 2.0)
        assert _select_ranges([r]) == [r]


class TestGetOptimalPdbsDf:

    def _make_df(self) -> pd.DataFrame:
        data = {
            "id": ["1ABC", "2XYZ", "3DEF"],
            "resolution": [2.0, 3.5, 1.5],
            "chains": ["A=1-180", "A=1-200", "A=1-100"],
            "type": ["PDB", "PDB", "PDB"],
        }
        df = pd.DataFrame(data)
        df[["chain_type", "chain_start", "chain_end"]] = df["chains"].str.extract(r"(.+)?=(\d+)-(\d+)")
        df["chain_type"] = df["chain_type"].apply(lambda x: str(x).split("/"))
        df["chain_start"] = df["chain_start"].astype(int)
        df["chain_end"] = df["chain_end"].astype(int)
        df["coverage"] = df["chain_end"] - df["chain_start"]
        df = df.set_index("id", drop=False)
        return df

    def test_returns_dataframe_subset(self):
        df = self._make_df()
        result = _get_optimal_pdbs_df(df, seq_length=200, min_res_val=2.5)
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        assert all(col in result.columns for col in ["id", "resolution", "coverage"])

    def test_all_ids_are_valid(self):
        df = self._make_df()
        result = _get_optimal_pdbs_df(df, seq_length=200, min_res_val=2.5)
        valid_ids = {"1ABC", "2XYZ", "3DEF"}
        assert all(pdb_id in valid_ids for pdb_id in result["id"].tolist())


class TestPdbCleanPrimitive:
    ENTRY_ID = "P02144"
    GENE_NAME = "MB"
    SCAN_ID = "integration_scan_001"
    OUTPUT = "integration_output"

    def test_pdb_clean_returns_datastore_address(self, disk_datastore, tmp_path, monkeypatch):
        """pdb_clean returns a valid datastore address pointing to a JSON summary."""
        monkeypatch.setenv("PROTEOMESCAN_CACHE_ROOT", str(tmp_path / "ps_cache"))
        config.set_datastore(disk_datastore)
        result_address = pdb_clean(
            gene_name=self.GENE_NAME,
            entry_id=self.ENTRY_ID,
            scan_id=self.SCAN_ID,
            output=self.OUTPUT,
            min_res_val=2.5,
        )

        assert result_address is not None
        assert result_address.startswith("deepchem://")

    def test_pdb_clean_summary_structure(self, disk_datastore, tmp_path, monkeypatch):
        """JSON summary contains required keys with correct types."""
        monkeypatch.setenv("PROTEOMESCAN_CACHE_ROOT", str(tmp_path / "ps_cache"))
        config.set_datastore(disk_datastore)
        result_address = pdb_clean(
            gene_name=self.GENE_NAME,
            entry_id=self.ENTRY_ID,
            scan_id=self.SCAN_ID,
            output=self.OUTPUT,
        )

        summary_raw = disk_datastore.get(result_address)
        summary = json.loads(summary_raw) if isinstance(summary_raw, str) else summary_raw

        assert "gene_name" in summary
        assert summary["gene_name"] == self.GENE_NAME
        assert "entry_id" in summary
        assert summary["entry_id"] == self.ENTRY_ID
        assert "scan_id" in summary
        assert summary["scan_id"] == self.SCAN_ID
        assert "selected_pdb_ids" in summary
        assert isinstance(summary["selected_pdb_ids"], list)
        assert len(summary["selected_pdb_ids"]) > 0
        assert "cleaned_pdb_addresses" in summary
        assert isinstance(summary["cleaned_pdb_addresses"], dict)
        assert "pdbs_metadata_csv_address" in summary
        assert isinstance(summary["pdbs_metadata_csv_address"], str)
        assert summary["pdbs_metadata_csv_address"].startswith("deepchem://")

    def test_cleaned_pdb_addresses_exist_in_datastore(self, disk_datastore, tmp_path, monkeypatch):
        """Every address in cleaned_pdb_addresses resolves to a real datastore object."""
        monkeypatch.setenv("PROTEOMESCAN_CACHE_ROOT", str(tmp_path / "ps_cache"))
        config.set_datastore(disk_datastore)
        result_address = pdb_clean(
            gene_name=self.GENE_NAME,
            entry_id=self.ENTRY_ID,
            scan_id=self.SCAN_ID,
            output=self.OUTPUT,
        )

        summary_raw = disk_datastore.get(result_address)
        summary = json.loads(summary_raw) if isinstance(summary_raw, str) else summary_raw

        for pdb_id, addr in summary["cleaned_pdb_addresses"].items():
            assert addr.startswith("deepchem://"), f"{pdb_id} address is not a deepchem:// address"
            assert disk_datastore.exists(addr), f"datastore object missing for {pdb_id}: {addr}"

    def test_pdb_clean_raises_without_datastore(self, monkeypatch, tmp_path):
        """pdb_clean raises ValueError when no datastore is configured."""
        monkeypatch.setenv("PROTEOMESCAN_CACHE_ROOT", str(tmp_path / "ps_cache"))
        config.set_datastore(None)

        with pytest.raises(ValueError, match="Datastore not configured"):
            pdb_clean(
                gene_name=self.GENE_NAME,
                entry_id=self.ENTRY_ID,
                scan_id=self.SCAN_ID,
                output=self.OUTPUT,
            )

    def test_pdb_cleaner_produces_output(self, tmp_path, monkeypatch):
        """_pdb_cleaner writes a non-empty cleaned PDB file."""
        from deepchem_server.core.primitives.proteome_scan import cache as ps_cache
        from deepchem_server.core.primitives.proteome_scan.pdb_clean import _pdb_cleaner

        scan_id = "integ_test"
        gene_name = "TESTGENE"
        pdb_id = "1ABC"
        monkeypatch.setenv("PROTEOMESCAN_CACHE_ROOT", str(tmp_path / "ps_cache"))

        raw_path = ps_cache.pdb_raw_path(scan_id, gene_name, pdb_id)
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(str(MINIMAL_PDB), str(raw_path))

        result = _pdb_cleaner(gene_name, pdb_id, scan_id, remove_chains=[])
        assert result is not None
        assert result.exists()
        assert result.stat().st_size > 0
