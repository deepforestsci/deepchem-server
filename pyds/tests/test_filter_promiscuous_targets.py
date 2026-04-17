"""Tests for FilterPromiscuousTargets primitive."""

import json
import os
import tempfile
import time
import uuid

import pandas as pd
import pytest
import responses

from pyds.data import Data
from pyds.primitives import FilterPromiscuousTargets
from pyds.settings import Settings


class TestFilterPromiscuousTargetsUnit:
    """Unit tests with mocked HTTP responses."""

    @responses.activate
    def test_init(self, test_settings: Settings) -> None:
        """Should initialize with settings."""
        client = FilterPromiscuousTargets(settings=test_settings)
        assert client.settings == test_settings
        assert client.base_url == "http://localhost:8000"

    @responses.activate
    def test_run_success(self, test_settings: Settings) -> None:
        """Should return filter_results_address on success."""
        mock_response = {"filter_results_address": "deepchem://test_profile/test_project/results.json"}
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/filter-promiscuous-targets",
            json=mock_response,
            status=200,
        )

        client = FilterPromiscuousTargets(settings=test_settings)
        result = client.run(
            scan_result_addresses=[
                "deepchem://test_profile/test_project/scan_a.csv",
                "deepchem://test_profile/test_project/scan_b.csv",
            ],
            thresholds=[[15, 1]],
            output="test_output",
        )

        assert result == mock_response

    @responses.activate
    def test_run_multiple_thresholds(self, test_settings: Settings) -> None:
        """Should support multiple threshold configurations."""
        mock_response = {"filter_results_address": "deepchem://test_profile/test_project/results.json"}
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/filter-promiscuous-targets",
            json=mock_response,
            status=200,
        )

        client = FilterPromiscuousTargets(settings=test_settings)
        client.run(
            scan_result_addresses=["deepchem://test_profile/test_project/scan.csv"],
            thresholds=[[15, 1], [25, 2]],
            output="test_output",
        )

        body = json.loads(responses.calls[0].request.body or "{}")
        assert body["thresholds"] == [[15, 1], [25, 2]]

    def test_run_missing_settings(self, test_settings_not_configured: Settings) -> None:
        """Should raise when settings are incomplete."""
        client = FilterPromiscuousTargets(settings=test_settings_not_configured)
        with pytest.raises(ValueError, match="Missing required settings"):
            client.run(
                scan_result_addresses=["deepchem://addr/scan.csv"],
                thresholds=[[15, 1]],
                output="test_output",
            )

    @responses.activate
    def test_run_server_error(self, test_settings: Settings) -> None:
        """Should raise on server error."""
        responses.add(
            responses.POST,
            "http://localhost:8000/primitive/filter-promiscuous-targets",
            json={"detail": "filter failed"},
            status=500,
        )

        client = FilterPromiscuousTargets(settings=test_settings)
        with pytest.raises(Exception):
            client.run(
                scan_result_addresses=["deepchem://addr/scan.csv"],
                thresholds=[[15, 1]],
                output="test_output",
            )


class TestFilterPromiscuousTargetsLive:
    """Integration tests against a live server at localhost:8000."""

    @pytest.fixture
    def sample_data_a(self) -> dict[str, list]:
        """Sample scan data with GENE_A and GENE_B in top 50%."""
        return {
            "gene_name": ["GENE_A", "GENE_B", "GENE_C", "GENE_D"],
            "top_score": [-9.0, -8.5, -7.0, -6.0],
        }

    @pytest.fixture
    def sample_data_b(self) -> dict[str, list]:
        """Sample scan data with GENE_A and GENE_E in top 50%."""
        return {
            "gene_name": ["GENE_A", "GENE_E", "GENE_F", "GENE_G"],
            "top_score": [-8.8, -8.0, -7.5, -6.5],
        }

    def _unique_name(self, prefix: str) -> str:
        """Generate unique output prefix."""
        return f"{prefix}_{uuid.uuid4().hex[:8]}_{int(time.time())}"

    def _upload_csv(self, data: dict[str, list], filename: str, client: Data) -> str:
        """Upload data as CSV and return dataset address."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            tmp_path = f.name

        try:
            pd.DataFrame(data).to_csv(tmp_path, index=False)
            result = client.upload_data(file_path=tmp_path, filename=filename)
            return result["dataset_address"]
        finally:
            os.unlink(tmp_path)

    def _download_json(self, address: str, client: Data) -> dict:
        """Download and parse JSON from datastore."""
        filename = address.split("/")[-1]
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            tmp_path = f.name

        try:
            client.download_data(address=filename, destination_path=tmp_path)
            with open(tmp_path) as f:
                return json.load(f)
        finally:
            os.unlink(tmp_path)

    def _download_csv(self, address: str, client: Data) -> pd.DataFrame:
        """Download and parse CSV from datastore."""
        filename = address.split("/")[-1]
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            tmp_path = f.name

        try:
            client.download_data(address=filename, destination_path=tmp_path)
            return pd.read_csv(tmp_path)
        finally:
            os.unlink(tmp_path)

    def test_filter_basic(
        self,
        live_filter_promiscuous_targets_client: FilterPromiscuousTargets,
        live_data_client: Data,
        sample_data_a: dict[str, list],
        sample_data_b: dict[str, list],
    ) -> None:
        """Should return valid filter_results_address."""
        uid = uuid.uuid4().hex[:8]
        addr_a = self._upload_csv(sample_data_a, f"scan_a_{uid}.csv", live_data_client)
        addr_b = self._upload_csv(sample_data_b, f"scan_b_{uid}.csv", live_data_client)

        result = live_filter_promiscuous_targets_client.run(
            scan_result_addresses=[addr_a, addr_b],
            thresholds=[[50, 2]],
            output=self._unique_name("test"),
        )

        assert "filter_results_address" in result
        assert result["filter_results_address"].endswith("_results.json")

    def test_promiscuous_gene_identified(
        self,
        live_filter_promiscuous_targets_client: FilterPromiscuousTargets,
        live_data_client: Data,
        sample_data_a: dict[str, list],
        sample_data_b: dict[str, list],
    ) -> None:
        """Should identify GENE_A as promiscuous at threshold (50, 2)."""
        uid = uuid.uuid4().hex[:8]
        addr_a = self._upload_csv(sample_data_a, f"scan_a_{uid}.csv", live_data_client)
        addr_b = self._upload_csv(sample_data_b, f"scan_b_{uid}.csv", live_data_client)

        result = live_filter_promiscuous_targets_client.run(
            scan_result_addresses=[addr_a, addr_b],
            thresholds=[[50, 2]],
            output=self._unique_name("test"),
        )

        results = self._download_json(result["filter_results_address"], live_data_client)
        promiscuous = results["promiscuous_targets"]["50%_2"]

        assert "GENE_A" in promiscuous
        assert "GENE_B" not in promiscuous
        assert "GENE_E" not in promiscuous

    def test_filtered_csv_content(
        self,
        live_filter_promiscuous_targets_client: FilterPromiscuousTargets,
        live_data_client: Data,
        sample_data_a: dict[str, list],
        sample_data_b: dict[str, list],
    ) -> None:
        """Should exclude promiscuous genes from filtered CSVs."""
        uid = uuid.uuid4().hex[:8]
        addr_a = self._upload_csv(sample_data_a, f"scan_a_{uid}.csv", live_data_client)
        addr_b = self._upload_csv(sample_data_b, f"scan_b_{uid}.csv", live_data_client)

        result = live_filter_promiscuous_targets_client.run(
            scan_result_addresses=[addr_a, addr_b],
            thresholds=[[50, 2]],
            output=self._unique_name("test"),
        )

        results = self._download_json(result["filter_results_address"], live_data_client)

        for filtered_addr in results["filtered_results"]["50%_2"].values():
            df = self._download_csv(filtered_addr, live_data_client)
            assert "GENE_A" not in df["gene_name"].values

    def test_multiple_thresholds(
        self,
        live_filter_promiscuous_targets_client: FilterPromiscuousTargets,
        live_data_client: Data,
        sample_data_a: dict[str, list],
        sample_data_b: dict[str, list],
    ) -> None:
        """Should handle multiple thresholds independently."""
        uid = uuid.uuid4().hex[:8]
        addr_a = self._upload_csv(sample_data_a, f"scan_a_{uid}.csv", live_data_client)
        addr_b = self._upload_csv(sample_data_b, f"scan_b_{uid}.csv", live_data_client)

        result = live_filter_promiscuous_targets_client.run(
            scan_result_addresses=[addr_a, addr_b],
            thresholds=[[50, 2], [50, 1]],
            output=self._unique_name("test"),
        )

        results = self._download_json(result["filter_results_address"], live_data_client)

        assert "50%_2" in results["filtered_results"]
        assert "50%_1" in results["filtered_results"]
        assert "50%_2" in results["promiscuous_targets"]
        assert "50%_1" in results["promiscuous_targets"]

    def test_missing_gene_name_column(
        self,
        live_filter_promiscuous_targets_client: FilterPromiscuousTargets,
        live_data_client: Data,
    ) -> None:
        """Should raise when CSV is missing required column."""
        uid = uuid.uuid4().hex[:8]
        bad_data = {"protein": ["P1", "P2"], "score": [1.0, 2.0]}
        addr_bad = self._upload_csv(bad_data, f"bad_{uid}.csv", live_data_client)  # type: ignore

        with pytest.raises(Exception):
            live_filter_promiscuous_targets_client.run(
                scan_result_addresses=[addr_bad],
                thresholds=[[50, 1]],
                output=self._unique_name("test"),
            )
