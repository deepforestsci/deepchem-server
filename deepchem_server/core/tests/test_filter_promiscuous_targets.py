import pandas as pd
import pytest

from deepchem_server.core import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.primitives.proteome_scan.filter_promiscuous_targets import (
    filter_promiscuous_targets,
)


class TestFilterPromiscuousTargets:
    """Tests for the filter_promiscuous_targets primitive."""

    @pytest.fixture
    def sample_data_a(self) -> list[tuple[str, float]]:
        """Sample data for test case A - 4 rows with GENE_A in top 2."""
        return [
            ("GENE_A", -9.0),
            ("GENE_B", -8.5),
            ("GENE_C", -7.0),
            ("GENE_D", -6.0),
        ]

    @pytest.fixture
    def sample_data_b(self) -> list[tuple[str, float]]:
        """Sample data for test case B - 4 rows with GENE_A in top 2."""
        return [
            ("GENE_A", -8.8),
            ("GENE_E", -8.0),
            ("GENE_F", -7.5),
            ("GENE_G", -6.5),
        ]

    def _create_scan_csv(self, rows: list[tuple[str, float]], datastore, filename: str) -> str:
        """Helper to create a scan result CSV in the datastore."""
        df = pd.DataFrame({"gene_name": [r[0] for r in rows], "top_score": [r[1] for r in rows]})
        card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
        return datastore.upload_data_from_memory(df, filename, card)

    def test_basic_filtering(self, disk_datastore, sample_data_a, sample_data_b):
        """Should identify targets appearing in top-N of multiple inputs."""
        config.set_datastore(disk_datastore)

        addr_a = self._create_scan_csv(sample_data_a, disk_datastore, "scan_a.csv")
        addr_b = self._create_scan_csv(sample_data_b, disk_datastore, "scan_b.csv")

        result_addr = filter_promiscuous_targets(
            scan_result_addresses=[addr_a, addr_b],
            thresholds=[[50, 2]],
            output="test_output",
        )

        assert result_addr.startswith("deepchem://")
        assert result_addr.endswith("_results.json")

        results = disk_datastore.get(result_addr)
        assert "promiscuous_targets" in results
        assert "promiscuous_targets_address" in results
        assert "filtered_results" in results

        promiscuous = results["promiscuous_targets"]["50%_2"]
        assert "GENE_A" in promiscuous
        assert "GENE_B" not in promiscuous
        assert "GENE_E" not in promiscuous

    def test_filtered_output_excludes_promiscuous(self, disk_datastore, sample_data_a, sample_data_b):
        """Filtered outputs should not contain promiscuous targets."""
        config.set_datastore(disk_datastore)

        addr_a = self._create_scan_csv(sample_data_a, disk_datastore, "scan_a.csv")
        addr_b = self._create_scan_csv(sample_data_b, disk_datastore, "scan_b.csv")

        result_addr = filter_promiscuous_targets(
            scan_result_addresses=[addr_a, addr_b],
            thresholds=[[50, 2]],
            output="test_output",
        )

        results = disk_datastore.get(result_addr)
        filtered = results["filtered_results"]["50%_2"]

        for _, filtered_addr in filtered.items():
            df = disk_datastore.get(filtered_addr)
            assert isinstance(df, pd.DataFrame)
            assert "GENE_A" not in df["gene_name"].values

    def test_multiple_thresholds(self, disk_datastore, sample_data_a, sample_data_b):
        """Should handle multiple threshold configurations independently."""
        config.set_datastore(disk_datastore)

        addr_a = self._create_scan_csv(sample_data_a, disk_datastore, "scan_a.csv")
        addr_b = self._create_scan_csv(sample_data_b, disk_datastore, "scan_b.csv")

        result_addr = filter_promiscuous_targets(
            scan_result_addresses=[addr_a, addr_b],
            thresholds=[[50, 2], [50, 1]],
            output="test_output",
        )

        results = disk_datastore.get(result_addr)

        assert "50%_2" in results["promiscuous_targets"]
        assert "50%_1" in results["promiscuous_targets"]

        promiscuous_1 = results["promiscuous_targets"]["50%_1"]
        for gene in ("GENE_A", "GENE_B", "GENE_E"):
            assert gene in promiscuous_1

        assert "50%_2" in results["filtered_results"]
        assert "50%_1" in results["filtered_results"]

    def test_no_promiscuous_targets(self, disk_datastore):
        """Should handle case where no targets meet promiscuity criteria."""
        config.set_datastore(disk_datastore)

        data_x = [("GENE_X1", -9.0), ("GENE_X2", -8.0), ("GENE_X3", -7.0), ("GENE_X4", -6.0)]
        data_y = [("GENE_Y1", -9.0), ("GENE_Y2", -8.0), ("GENE_Y3", -7.0), ("GENE_Y4", -6.0)]

        addr_x = self._create_scan_csv(data_x, disk_datastore, "scan_x.csv")
        addr_y = self._create_scan_csv(data_y, disk_datastore, "scan_y.csv")

        result_addr = filter_promiscuous_targets(
            scan_result_addresses=[addr_x, addr_y],
            thresholds=[[50, 2]],
            output="test_output",
        )

        results = disk_datastore.get(result_addr)
        assert results["promiscuous_targets"]["50%_2"] == []

        for filtered_addr in results["filtered_results"]["50%_2"].values():
            df = disk_datastore.get(filtered_addr)
            assert len(df) == 4

    def test_missing_gene_name_column_raises(self, disk_datastore):
        """Should raise KeyError when required column is missing."""
        config.set_datastore(disk_datastore)

        bad_df = pd.DataFrame({"protein": ["val1", "val2"], "score": [1.0, 2.0]})
        card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
        bad_addr = disk_datastore.upload_data_from_memory(bad_df, "bad_scan.csv", card)

        with pytest.raises(KeyError, match="gene_name"):
            filter_promiscuous_targets(
                scan_result_addresses=[bad_addr],
                thresholds=[[50, 1]],
                output="test_output",
            )

    def test_empty_addresses_raises(self, disk_datastore):
        """Should raise ValueError when no scan addresses provided."""
        config.set_datastore(disk_datastore)

        with pytest.raises(ValueError, match="At least one scan result address"):
            filter_promiscuous_targets(
                scan_result_addresses=[],
                thresholds=[[50, 1]],
                output="test_output",
            )
