"""
Tests for the ``proteome_scan.pdb_clean`` primitive.

These tests focus on the deepchem-server boundary (server cache
layout, resumability, datastore uploads, and helper-function
behaviour). They do not hit UniProt / PDBe-KB or PDBFixer - those
codepaths are already covered by ProteomeScan's own tests and are
intentionally left to integration suites.
"""

import json
import os

import pandas as pd
import pytest

from deepchem_server.core import config
from deepchem_server.core.primitives.proteome_scan import cache as ps_cache
from deepchem_server.core.primitives.proteome_scan.pdb_clean import (
    Range,
    pdb_clean,
    select_ranges,
)


@pytest.fixture
def proteome_scan_cache(tmp_path, monkeypatch):
    """Redirect the proteome_scan cache root to a pytest tmp_path."""
    cache_root = tmp_path / "proteome_scan_cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    monkeypatch.setenv("PROTEOMESCAN_CACHE_ROOT", str(cache_root))
    return cache_root


def _write_dummy_pdb(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("HEADER    DUMMY PDB\n")
        f.write("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n")
        f.write("END\n")


class TestRangeHelpers:
    """Unit tests for the ``Range`` and ``select_ranges`` helpers."""

    def test_range_repr_and_coverage(self):
        r = Range("1ABC", 10, 50, 2.0)
        assert r.id == "1ABC"
        assert r.coverage == 40
        assert "1ABC" in repr(r)

    def test_select_ranges_non_overlapping(self):
        ranges = [
            Range("A", 1, 100, 2.0),
            Range("B", 150, 250, 1.5),
            Range("C", 300, 400, 2.5),
        ]
        selected = select_ranges(ranges)
        ids = sorted(r.id for r in selected)
        assert ids == ["A", "B", "C"]

    def test_select_ranges_overlap_prefers_better_coverage(self):
        ranges = [
            Range("A", 1, 50, 2.0),
            Range("B", 10, 80, 1.9),
            Range("C", 200, 300, 2.0),
        ]
        selected = select_ranges(ranges)
        ids = {r.id for r in selected}
        assert "C" in ids
        assert "A" in ids or "B" in ids


class TestPdbCleanResumability:
    """Tests that exercise the on-disk cache / resumability path."""

    def test_resumes_from_existing_csv(self, disk_datastore, proteome_scan_cache):
        """If the per-gene CSV exists on disk, the primitive must
        reuse it and only upload datastore artifacts."""
        config.set_datastore(disk_datastore)

        scan_id = "scan_resume"
        gene_name = "GENE1"
        entry_id = "P12345"

        gene_dir = ps_cache.get_gene_dir(scan_id, gene_name)
        cleaned_path = str(gene_dir / f"cleaned_g_{gene_name}_p_1ABC.pdb")
        _write_dummy_pdb(cleaned_path)

        df = pd.DataFrame(
            {
                "id": ["1ABC"],
                "chain_type": ["[A]"],
                "chain_start": [1],
                "chain_end": [100],
                "resolution": [1.5],
                "chains": ["A=1-100"],
                "coverage": [99],
                "path": [cleaned_path],
            }
        ).set_index("id", drop=False)
        df.to_csv(str(gene_dir / f"{gene_name}_pdbs.csv"))

        result_address = pdb_clean(
            gene_name=gene_name,
            entry_id=entry_id,
            scan_id=scan_id,
            output="resume_test",
            min_res_val=3.0,
        )

        assert result_address.startswith("deepchem://")
        assert result_address.endswith("_metadata.json")

        metadata = disk_datastore.get(result_address)
        if isinstance(metadata, str):
            metadata = json.loads(metadata)

        assert metadata["gene_name"] == gene_name
        assert metadata["scan_id"] == scan_id
        assert metadata["num_pdbs"] == 1
        assert metadata["selected_pdb_ids"] == ["1ABC"]

        csv_address = metadata["pdbs_csv_address"]
        assert csv_address.startswith("deepchem://")
        reloaded = disk_datastore.get(csv_address)
        if isinstance(reloaded, str):
            from io import StringIO

            reloaded = pd.read_csv(StringIO(reloaded))
        assert "path" in reloaded.columns
        assert len(reloaded) == 1

        pdb_addresses = metadata["pdb_addresses"]
        assert "1ABC" in pdb_addresses
        assert pdb_addresses["1ABC"].startswith("deepchem://")

    def test_cache_layout_under_scan_dir(self, proteome_scan_cache):
        """Cache helpers lay out files exactly like ProteomeScan."""
        scan_id = "scan_layout"
        gene = "GENE2"
        ligand = "LIG"
        pdb = "1XYZ"
        assert str(ps_cache.get_scan_dir(scan_id)).endswith(
            os.path.join("proteome_scan_cache", scan_id)
        )
        assert str(ps_cache.get_gene_dir(scan_id, gene)).endswith(os.path.join(scan_id, gene))
        assert str(ps_cache.get_ligand_dir(scan_id, ligand)).endswith(os.path.join(scan_id, ligand))
        assert str(ps_cache.raw_pdb_path(scan_id, gene, pdb)).endswith(f"g_{gene}_p_{pdb}.pdb")
        assert str(ps_cache.cleaned_pdb_path(scan_id, gene, pdb)).endswith(
            f"cleaned_g_{gene}_p_{pdb}.pdb"
        )
        assert str(ps_cache.complex_path(scan_id, gene, pdb, ligand)).endswith(
            os.path.join("complexes", f"complex_{gene}_{pdb}_{ligand}.pdb")
        )

    def test_missing_required_params(self, disk_datastore, proteome_scan_cache):
        config.set_datastore(disk_datastore)
        with pytest.raises(ValueError):
            pdb_clean(
                gene_name="",
                entry_id="P12345",
                scan_id="s",
                output="o",
            )
        with pytest.raises(ValueError):
            pdb_clean(
                gene_name="g",
                entry_id="",
                scan_id="s",
                output="o",
            )
        with pytest.raises(ValueError):
            pdb_clean(
                gene_name="g",
                entry_id="p",
                scan_id="",
                output="o",
            )
