"""
Unit tests for DelDenoise primitive using live server.
"""

import os
import time
import uuid

import pytest

from pyds.data import Data
from pyds.primitives import DelDenoise
from pyds.settings import Settings


ASSETS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")
DEL_CSV = os.path.join(ASSETS_DIR, "del_test_dataset.csv")


class TestDelDenoise:
    """Unit tests for DelDenoise primitive."""

    def test_init(self, test_settings: Settings) -> None:
        """Test DelDenoise initialization."""
        client = DelDenoise(settings=test_settings)
        assert client.settings == test_settings
        assert client.base_url == "http://localhost:8000"

    def test_run_unified(
        self,
        live_del_denoise_client: DelDenoise,
        live_data_client: Data,
    ) -> None:
        """Test unified strategy run on live server."""
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        upload_result = live_data_client.upload_data(
            file_path=DEL_CSV,
            filename=f"del_raw_{test_id}_{timestamp}.csv",
            description="DEL test data for unified denoising",
        )
        dataset_address = upload_result["dataset_address"]

        result = live_del_denoise_client.run(
            dataset_address=dataset_address,
            output_key=f"del_unified_{test_id}_{timestamp}",
            strategy="unified",
        )

        assert "denoised_dataset_address" in result

    def test_run_unified_with_hit_labels(
        self,
        live_del_denoise_client: DelDenoise,
        live_data_client: Data,
    ) -> None:
        """Test unified strategy with hit labels on live server."""
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        upload_result = live_data_client.upload_data(
            file_path=DEL_CSV,
            filename=f"del_raw_hits_{test_id}_{timestamp}.csv",
            description="DEL test data for unified denoising with hits",
        )
        dataset_address = upload_result["dataset_address"]

        result = live_del_denoise_client.run(
            dataset_address=dataset_address,
            output_key=f"del_unified_hits_{test_id}_{timestamp}",
            strategy="unified",
            add_hit_labels=True,
            hit_percentile=90.0,
        )

        assert "denoised_dataset_address" in result

    def test_run_non_unified(
        self,
        live_del_denoise_client: DelDenoise,
        live_data_client: Data,
    ) -> None:
        """Test non-unified (z-score) strategy on live server."""
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        upload_result = live_data_client.upload_data(
            file_path=DEL_CSV,
            filename=f"del_raw_nu_{test_id}_{timestamp}.csv",
            description="DEL test data for non-unified denoising",
        )
        dataset_address = upload_result["dataset_address"]

        result = live_del_denoise_client.run(
            dataset_address=dataset_address,
            output_key=f"del_non_unified_{test_id}_{timestamp}",
            strategy="non_unified",
            add_hit_labels=True,
            hit_percentile=90.0,
        )

        assert "denoised_dataset_address" in result

    def test_run_with_disynthon_pairs(
        self,
        live_del_denoise_client: DelDenoise,
        live_data_client: Data,
    ) -> None:
        """Test disynthon collapse mode on live server."""
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        upload_result = live_data_client.upload_data(
            file_path=DEL_CSV,
            filename=f"del_raw_disynthon_{test_id}_{timestamp}.csv",
            description="DEL test data for disynthon denoising",
        )
        dataset_address = upload_result["dataset_address"]

        result = live_del_denoise_client.run(
            dataset_address=dataset_address,
            output_key=f"del_disynthon_{test_id}_{timestamp}",
            strategy="unified",
            use_disynthon_pairs=True,
        )

        assert "denoised_dataset_address" in result

    def test_run_with_custom_count_cols(
        self,
        live_del_denoise_client: DelDenoise,
        live_data_client: Data,
    ) -> None:
        """Test run with explicit control/target column names on live server."""
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        upload_result = live_data_client.upload_data(
            file_path=DEL_CSV,
            filename=f"del_raw_cols_{test_id}_{timestamp}.csv",
            description="DEL test data with custom columns",
        )
        dataset_address = upload_result["dataset_address"]

        result = live_del_denoise_client.run(
            dataset_address=dataset_address,
            output_key=f"del_custom_cols_{test_id}_{timestamp}",
            strategy="unified",
            control_cols=["seq_matrix_1", "seq_matrix_2", "seq_matrix_3"],
            target_cols=["seq_target_1", "seq_target_2", "seq_target_3"],
        )

        assert "denoised_dataset_address" in result

    def test_run_with_profile_project_override(
        self,
        live_del_denoise_client: DelDenoise,
        live_data_client: Data,
    ) -> None:
        """Test run with explicit profile and project override on live server."""
        test_id = str(uuid.uuid4())[:8]
        timestamp = str(int(time.time()))

        upload_result = live_data_client.upload_data(
            file_path=DEL_CSV,
            filename=f"del_raw_override_{test_id}_{timestamp}.csv",
            description="DEL test data with profile/project override",
        )
        dataset_address = upload_result["dataset_address"]

        result = live_del_denoise_client.run(
            dataset_address=dataset_address,
            output_key=f"del_override_{test_id}_{timestamp}",
            strategy="unified",
            profile_name="test_profile",
            project_name="test_project",
        )

        assert "denoised_dataset_address" in result

    def test_run_missing_settings(self, test_settings_not_configured: Settings) -> None:
        """Test run with missing profile/project settings."""
        client = DelDenoise(settings=test_settings_not_configured)
        with pytest.raises(ValueError, match="Missing required settings"):
            client.run(
                dataset_address="test/raw_del.csv",
                output_key="denoised",
            )
