"""
Tests for the ``proteome_scan.run_docking`` primitive.

These tests exercise the deepchem-server boundary and the
ProteomeScan-style on-disk layout. They stub out the VINA docking
call so the suite does not require a functional ``VinaPoseGenerator``
installation.
"""

import json
import os

import pandas as pd
import pytest

from deepchem_server.core import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.primitives.proteome_scan import cache as ps_cache
from deepchem_server.core.primitives.proteome_scan import docking as dk


@pytest.fixture
def proteome_scan_cache(tmp_path, monkeypatch):
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


@pytest.fixture
def ligand_address(disk_datastore):
    sdf = ("LIG\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
           "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
           "M  END\n$$$$\n")
    card = DataCard(address="", file_type="sdf", data_type="sdf")
    return disk_datastore.upload_data_from_memory(sdf, "test_lig.sdf", card)


@pytest.fixture
def prepared_gene(proteome_scan_cache):
    """Pre-populate <cache>/<scan>/<gene>/ with a cleaned PDB + CSV."""
    scan_id = "scan_dock"
    gene_name = "GENE1"
    pdb_id = "1ABC"
    gene_dir = ps_cache.get_gene_dir(scan_id, gene_name)
    cleaned_path = str(gene_dir / f"cleaned_g_{gene_name}_p_{pdb_id}.pdb")
    _write_dummy_pdb(cleaned_path)

    df = pd.DataFrame({
        "id": [pdb_id],
        "chain_type": ["[A]"],
        "chain_start": [1],
        "chain_end": [100],
        "resolution": [1.5],
        "chains": ["A=1-100"],
        "coverage": [99],
        "path": [cleaned_path],
    }).set_index("id", drop=False)
    df.to_csv(str(gene_dir / f"{gene_name}_pdbs.csv"))
    return {
        "scan_id": scan_id,
        "gene_name": gene_name,
        "pdb_id": pdb_id,
        "cleaned_path": cleaned_path,
    }


class TestRunDocking:

    def test_run_docking_writes_expected_layout(self, disk_datastore, ligand_address, prepared_gene, monkeypatch):
        """run_docking must write complexes + top_score CSV under the
        server cache and upload everything to the datastore."""
        config.set_datastore(disk_datastore)

        scan_id = prepared_gene["scan_id"]
        gene_name = prepared_gene["gene_name"]
        pdb_id = prepared_gene["pdb_id"]
        ligand_name = "LIG"

        def _fake_vina(pdb_path, lig_path, work_dir, exhaustiveness=32, num_modes=8):
            assert os.path.exists(pdb_path)
            assert os.path.exists(lig_path)
            return ("tuple",), [-7.5, -7.2]

        def _fake_write_complex(complex_tuple, output_path):
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, "w") as f:
                f.write("HEADER    FAKE COMPLEX\nEND\n")
            return True

        monkeypatch.setattr(dk, "_vina_docking", _fake_vina)
        monkeypatch.setattr(dk, "_write_complex_pdb", _fake_write_complex)

        result_addr = dk.run_docking(
            gene_name=gene_name,
            ligand_name=ligand_name,
            ligand_address=ligand_address,
            scan_id=scan_id,
            output="dock_run_1",
            exhaustiveness=1,
            num_modes=1,
        )

        assert result_addr.startswith("deepchem://")
        summary = disk_datastore.get(result_addr)
        if isinstance(summary, str):
            summary = json.loads(summary)

        assert summary["gene_name"] == gene_name
        assert summary["ligand_name"] == ligand_name
        assert summary["scan_id"] == scan_id
        assert summary["top_score_csv_address"].startswith("deepchem://")

        top_csv = ps_cache.top_score_gene_ligand_csv_path(scan_id, gene_name, ligand_name)
        assert top_csv.exists()
        csv_df = pd.read_csv(top_csv)
        assert "top_score" in csv_df.columns
        assert "gene_name" in csv_df.columns
        assert (csv_df["gene_name"] == gene_name).all()

        complex_path = ps_cache.complex_path(scan_id, gene_name, pdb_id, ligand_name)
        assert complex_path.exists()

        assert "complex_addresses" in summary
        assert any(complex_path.stem == k for k in summary["complex_addresses"])

    def test_run_docking_resumes_when_csv_exists(self, disk_datastore, ligand_address, prepared_gene, monkeypatch):
        config.set_datastore(disk_datastore)
        scan_id = prepared_gene["scan_id"]
        gene_name = prepared_gene["gene_name"]
        ligand_name = "LIG"

        top_csv = ps_cache.top_score_gene_ligand_csv_path(scan_id, gene_name, ligand_name)
        top_csv.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({
            "id": ["1ABC"],
            "chains": ["A=1-100"],
            "resolution": [1.5],
            "coverage": [99],
            "top_score": [-9.0],
            "scores": [[-9.0]],
            "gene_name": [gene_name],
        }).to_csv(top_csv, index=False)

        called = {"n": 0}

        def _should_not_run(*args, **kwargs):
            called["n"] += 1
            raise AssertionError("VINA docking must not run on resume path")

        monkeypatch.setattr(dk, "_vina_docking", _should_not_run)

        result_addr = dk.run_docking(
            gene_name=gene_name,
            ligand_name=ligand_name,
            ligand_address=ligand_address,
            scan_id=scan_id,
            output="dock_resume",
        )

        assert called["n"] == 0
        summary = disk_datastore.get(result_addr)
        if isinstance(summary, str):
            summary = json.loads(summary)
        assert summary["gene_name"] == gene_name
        assert summary["top_score_csv_address"].startswith("deepchem://")

    def test_run_docking_requires_gene_csv(self, disk_datastore, ligand_address, proteome_scan_cache):
        config.set_datastore(disk_datastore)
        with pytest.raises(FileNotFoundError):
            dk.run_docking(
                gene_name="NONE",
                ligand_name="LIG",
                ligand_address=ligand_address,
                scan_id="no_scan",
                output="dock_fail",
            )
