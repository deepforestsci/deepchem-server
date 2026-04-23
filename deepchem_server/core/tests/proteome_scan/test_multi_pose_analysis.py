"""
Tests for ``proteome_scan.run_multi_pose_analysis``.

These tests assume PyMOL is available and fpocket is runnable via Docker
(``fpocket/fpocket`` image). They validate both the deepchem-server
boundary (datastore/cache interactions) and the end-to-end pose analysis
path.
"""

import json
import os

import pandas as pd
import pytest

from deepchem_server.core import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.primitives.proteome_scan import cache
from deepchem_server.core.primitives.proteome_scan import multi_pose_analysis


@pytest.fixture
def proteome_scan_cache(tmp_path, monkeypatch):
    cache_root = tmp_path / "proteome_scan_cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    monkeypatch.setenv("PROTEOMESCAN_CACHE_ROOT", str(cache_root))
    return cache_root


@pytest.fixture
def sample_complex(disk_datastore):
    pdb_content = (
        "HEADER    COMPLEX\n"
        "ATOM      1  N   SER A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "ATOM      2  CA  SER A   1      11.000  10.000  10.000  1.00 20.00           C\n"
        "HETATM    3  C1  LIG A 100      15.000  15.000  15.000  1.00 20.00           C\n"
        "END\n"
    )
    card = DataCard(address="", file_type="pdb", data_type="text/plain")
    return disk_datastore.upload_data_from_memory(pdb_content, "test_complex.pdb", card)


class TestArgValidation:

    def test_empty_complex_addresses_raises(self, disk_datastore, proteome_scan_cache):
        config.set_datastore(disk_datastore)
        with pytest.raises(ValueError, match="complex_addresses"):
            multi_pose_analysis.run_multi_pose_analysis(complex_addresses=[], output="o")

    def test_missing_output_raises(self, disk_datastore, proteome_scan_cache):
        config.set_datastore(disk_datastore)
        with pytest.raises(ValueError, match="output"):
            multi_pose_analysis.run_multi_pose_analysis(complex_addresses=["a"], output="")

    def test_no_datastore_raises(self, proteome_scan_cache):
        config.set_datastore(None)
        with pytest.raises(ValueError, match="Datastore"):
            multi_pose_analysis.run_multi_pose_analysis(complex_addresses=["a"], output="o")


class TestParsePocketData:

    def test_parse_pocket_data(self, tmp_path):
        info = tmp_path / "protein_info.txt"
        info.write_text(
            "Pocket 1 :\n"
            "Druggability Score : 0.234\n"
            "Number of Alpha Spheres : 15\n"
            "Volume : 234.56\n"
            "Mean local hydrophobic density : 12.34\n"
            "\n"
            "Pocket 2 :\n"
            "Druggability Score : 0.123\n"
            "Number of Alpha Spheres : 10\n"
            "Volume : 100.0\n"
            "Mean local hydrophobic density : 5.0\n"
        )

        pockets = multi_pose_analysis._parse_pocket_data(str(info))
        assert len(pockets) == 2
        assert pockets[0]["pocket_id"] == 1
        assert pockets[1]["pocket_id"] == 2
        assert pockets[0]["Druggability Score"] == 0.234
        assert pockets[0]["Number of Alpha Spheres"] == 15


class TestAggregationHelpers:

    def test_overall_caps_at_100(self):
        df = pd.DataFrame({"% Ligand inside pocket": [60.0, 50.0, 20.0]})
        assert multi_pose_analysis._get_overall_ligand_interactions(df) == 100

    def test_overall_sums_when_below_100(self):
        df = pd.DataFrame({"% Ligand inside pocket": [10.0, 20.0]})
        assert multi_pose_analysis._get_overall_ligand_interactions(df) == 30.0

    def test_topn_bucket_percentages(self):
        df = pd.DataFrame(
            {
                "% Ligand inside pocket": [30.0, 20.0, 15.0, 5.0],
                "Pocket_Druggability_Rank": [1, 2, 3, 4],
            }
        )
        assert multi_pose_analysis._get_total_top_n_bucket_percentages(df, 1) == 30.0
        assert multi_pose_analysis._get_total_top_n_bucket_percentages(df, 5) == 70.0


class TestEndToEnd:

    def test_run_multi_pose_analysis_basic(
        self, disk_datastore, proteome_scan_cache, sample_complex
    ):
        config.set_datastore(disk_datastore)
        result_addr = multi_pose_analysis.run_multi_pose_analysis(
            complex_addresses=[sample_complex],
            output="test_pose_analysis",
            num_processes=1,
            is_clean_up=True,
            scan_id="scan_mpa",
        )
        assert result_addr.startswith("deepchem://")
        summary = disk_datastore.get(result_addr)
        if isinstance(summary, str):
            summary = json.loads(summary)
        assert summary["num_complexes"] == 1
        assert "results_csv_address" in summary

        csv_addr = summary["results_csv_address"]
        df = disk_datastore.get(csv_addr)
        if isinstance(df, str):
            from io import StringIO

            df = pd.read_csv(StringIO(df))
        assert "complex" in df.columns
        for n in [1, 5, 10]:
            assert f"total % Ligand inside pockets (top{n} pockets)" in df.columns


class TestCacheLayout:
    """Verify that pose-analysis downloads land in the expected cache
    sub-directory even when the analyse_pose invocation itself is not
    runnable in this environment."""

    def test_download_to_cache_uses_scan_dir(
        self, disk_datastore, proteome_scan_cache, sample_complex
    ):
        config.set_datastore(disk_datastore)
        local = multi_pose_analysis._download_complex_to_cache(
            sample_complex, scan_id="scan_cache", output="ignored"
        )
        assert os.path.exists(local)
        assert "scan_cache" in local
        assert local.endswith(".pdb")

    def test_download_to_cache_without_scan_id(
        self, disk_datastore, proteome_scan_cache, sample_complex
    ):
        config.set_datastore(disk_datastore)
        local = multi_pose_analysis._download_complex_to_cache(
            sample_complex, scan_id=None, output="no_scan"
        )
        assert os.path.exists(local)
        root = str(cache.get_cache_root())
        assert local.startswith(root)
