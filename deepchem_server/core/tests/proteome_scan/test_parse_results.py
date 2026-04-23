"""
Tests for the proteome_scan.parse_results primitive.
"""

import json
from io import StringIO

import pandas as pd
import pytest

from deepchem_server.core import config
from deepchem_server.core.primitives.proteome_scan import cache as ps_cache
from deepchem_server.core.primitives.proteome_scan.parse_results import (
    parse_results,
)


@pytest.fixture
def proteome_scan_cache(tmp_path, monkeypatch):
    cache_root = tmp_path / "proteome_scan_cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    monkeypatch.setenv("PROTEOMESCAN_CACHE_ROOT", str(cache_root))
    return cache_root


class TestParseResults:
    def test_parse_results_aggregates_and_sorts(self, disk_datastore, proteome_scan_cache):
        """parse_results must concat per-gene CSVs, sort by top_score
        and dedupe by gene_name."""
        config.set_datastore(disk_datastore)

        scan_id = "scan_parse"
        ligand = "LIG"
        ligand_dir = ps_cache.get_ligand_dir(scan_id, ligand)

        df1 = pd.DataFrame(
            {
                "": ["1ABC", "1XYZ"],
                "chains": ["A=1-100", "A=1-100"],
                "resolution": [1.5, 1.8],
                "coverage": [99, 99],
                "top_score": [-9.2, -7.5],
                "scores": [[-9.2], [-7.5]],
                "gene_name": ["GENE1", "GENE1"],
            }
        )
        df1.to_csv(str(ligand_dir / f"top_score_GENE1_{ligand}.csv"), index=False)

        df2 = pd.DataFrame(
            {
                "": ["2DEF"],
                "chains": ["A=1-200"],
                "resolution": [1.2],
                "coverage": [199],
                "top_score": [-8.5],
                "scores": [[-8.5]],
                "gene_name": ["GENE2"],
            }
        )
        df2.to_csv(str(ligand_dir / f"top_score_GENE2_{ligand}.csv"), index=False)

        summary_addr = parse_results(
            scan_id=scan_id,
            ligands=[ligand],
            output="parse_run",
        )

        assert summary_addr.startswith("deepchem://")
        summary = disk_datastore.get(summary_addr)
        if isinstance(summary, str):
            summary = json.loads(summary)

        assert summary["scan_id"] == scan_id
        assert ligand in summary["aggregated"]
        agg = summary["aggregated"][ligand]
        assert agg["num_rows"] == 2

        out_csv = ps_cache.top_score_ligand_csv_path(scan_id, ligand)
        assert out_csv.exists()
        out_df = pd.read_csv(out_csv)
        assert list(out_df["gene_name"]) == ["GENE1", "GENE2"]
        assert out_df["top_score"].iloc[0] == -9.2

        csv_addr = agg["address"]
        reloaded = disk_datastore.get(csv_addr)
        if isinstance(reloaded, str):
            reloaded = pd.read_csv(StringIO(reloaded))
        pd.testing.assert_frame_equal(
            out_df.reset_index(drop=True),
            reloaded.reset_index(drop=True),
            check_dtype=False,
        )

    def test_parse_results_empty_ligand_folder(self, disk_datastore, proteome_scan_cache):
        config.set_datastore(disk_datastore)
        scan_id = "scan_parse_empty"
        ps_cache.get_ligand_dir(scan_id, "EMPTY")

        summary_addr = parse_results(
            scan_id=scan_id,
            ligands=["EMPTY"],
            output="parse_empty",
        )
        summary = disk_datastore.get(summary_addr)
        if isinstance(summary, str):
            summary = json.loads(summary)
        assert summary["aggregated"]["EMPTY"]["num_rows"] == 0
        assert summary["aggregated"]["EMPTY"]["address"] is None

    def test_parse_results_requires_args(self, disk_datastore, proteome_scan_cache):
        config.set_datastore(disk_datastore)
        with pytest.raises(ValueError):
            parse_results(scan_id="", ligands=["LIG"], output="o")
        with pytest.raises(ValueError):
            parse_results(scan_id="s", ligands=[], output="o")
        with pytest.raises(ValueError):
            parse_results(scan_id="s", ligands=["L"], output="")
