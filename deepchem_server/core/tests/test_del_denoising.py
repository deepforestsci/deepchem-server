"""Tests for DEL denoising primitive"""

import os

import numpy as np
import pandas as pd

from deepchem_server.core.common import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.primitives.del_denoising import (
    _calculate_normalized_enrichment_score,
    _calculate_poisson_enrichment,
    _poissfit,
    del_denoise,
)


ASSETS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")
DEL_CSV = os.path.join(ASSETS_DIR, "del_test_dataset.csv")


def _upload_del_csv(disk_datastore):
    """Helper to upload the DEL test CSV and return its address."""
    df = pd.read_csv(DEL_CSV)
    card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    return disk_datastore.upload_data_from_memory(df, "del_raw.csv", card)


def test_poissfit_basic():
    """Test Poisson CI on a known vector."""
    vec = pd.Series([10, 20, 30])
    lower, upper = _poissfit(vec, alpha=0.05)
    assert lower > 0
    assert upper > lower
    assert abs(lower - 15.262) < 0.01
    assert abs(upper - 25.744) < 0.01


def test_normalized_enrichment_score():
    """Match the expected z-score from deepchem-del test suite."""
    df = pd.DataFrame(
        {"seq_target_1": [10, 20, 30], "seq_target_2": [40, 50, 60], "seq_target_3": [70, 80, 90]}
    )
    row = df.iloc[0]
    total_sum = 60
    row_count = 3
    result = _calculate_normalized_enrichment_score(row, total_sum, row_count, "seq_target_1")
    assert abs(result - (-0.35355339059327373)) < 1e-6


def test_poisson_enrichment_matches_reference():
    """Ensure Poisson enrichment matches the pre-computed target_enrichment column."""
    df = pd.read_csv(DEL_CSV)
    control_cols = ["seq_matrix_1", "seq_matrix_2", "seq_matrix_3"]
    target_cols = ["seq_target_1", "seq_target_2", "seq_target_3"]
    result = _calculate_poisson_enrichment(df, control_cols, target_cols)
    np.testing.assert_allclose(
        result["Poisson_Enrichment"].values, df["target_enrichment"].values, rtol=1e-5
    )


def test_del_denoise_unified_regression(disk_datastore):
    """Unified strategy produces Poisson_Enrichment column."""
    config.set_datastore(disk_datastore)
    address = _upload_del_csv(disk_datastore)

    result_address = del_denoise(
        dataset_address=address, output_key="denoised_unified", strategy="unified"
    )

    result_df = disk_datastore.get(result_address)
    assert "Poisson_Enrichment" in result_df.columns
    assert len(result_df) > 0
    assert result_df["Poisson_Enrichment"].notna().all()


def test_del_denoise_unified_matches_reference(disk_datastore):
    """Unified enrichment values match the reference target_enrichment column."""
    config.set_datastore(disk_datastore)
    address = _upload_del_csv(disk_datastore)

    result_address = del_denoise(
        dataset_address=address, output_key="denoised_unified_check", strategy="unified"
    )

    result_df = disk_datastore.get(result_address)
    ref_df = pd.read_csv(DEL_CSV)
    ref_df = ref_df.dropna(subset=["smiles"]).drop_duplicates(subset=["smiles"])

    np.testing.assert_allclose(
        result_df["Poisson_Enrichment"].values, ref_df["target_enrichment"].values, rtol=1e-5
    )


def test_del_denoise_unified_classification(disk_datastore):
    """Unified + hit labels produces binary hits column."""
    config.set_datastore(disk_datastore)
    address = _upload_del_csv(disk_datastore)

    result_address = del_denoise(
        dataset_address=address,
        output_key="denoised_unified_cls",
        strategy="unified",
        add_hit_labels=True,
        hit_percentile=90.0,
    )

    result_df = disk_datastore.get(result_address)
    assert "Poisson_Enrichment" in result_df.columns
    assert "hits" in result_df.columns
    assert set(result_df["hits"].unique()).issubset({0, 1})


def test_del_denoise_non_unified(disk_datastore):
    """Non-unified strategy produces target and control enrichment scores."""
    config.set_datastore(disk_datastore)
    address = _upload_del_csv(disk_datastore)

    result_address = del_denoise(
        dataset_address=address, output_key="denoised_non_unified", strategy="non_unified"
    )

    result_df = disk_datastore.get(result_address)
    assert "Target_Enrichment_Score" in result_df.columns
    assert "Control_Enrichment_Score" in result_df.columns
    assert "seq_target_sum" in result_df.columns
    assert "seq_control_sum" in result_df.columns


def test_del_denoise_non_unified_with_hits(disk_datastore):
    """Non-unified + hit labels produces target_hits and control_hits."""
    config.set_datastore(disk_datastore)
    address = _upload_del_csv(disk_datastore)

    result_address = del_denoise(
        dataset_address=address,
        output_key="denoised_non_unified_cls",
        strategy="non_unified",
        add_hit_labels=True,
        hit_percentile=90.0,
    )

    result_df = disk_datastore.get(result_address)
    assert "target_hits" in result_df.columns
    assert "control_hits" in result_df.columns
    assert set(result_df["target_hits"].unique()).issubset({0, 1})
    assert set(result_df["control_hits"].unique()).issubset({0, 1})


def test_del_denoise_card_metadata(disk_datastore):
    """DataCard contains expected DEL metadata fields."""
    config.set_datastore(disk_datastore)
    address = _upload_del_csv(disk_datastore)

    result_address = del_denoise(
        dataset_address=address, output_key="denoised_card_test", strategy="unified"
    )

    card = disk_datastore.get_card(result_address, kind="data")
    assert card.strategy == "unified"
    assert card.parent == address
    assert card.n_input_rows > 0
    assert card.n_output_rows > 0


def test_del_denoise_with_disynthon_unified(disk_datastore):
    """use_disynthon_pairs=True with unified strategy produces disynthons + Poisson enrichment."""
    config.set_datastore(disk_datastore)
    address = _upload_del_csv(disk_datastore)

    result_address = del_denoise(
        dataset_address=address,
        output_key="denoised_disynthon_unified",
        strategy="unified",
        use_disynthon_pairs=True,
    )

    result_df = disk_datastore.get(result_address)
    assert "disynthons" in result_df.columns
    assert "Poisson_Enrichment" in result_df.columns
    assert len(result_df) > 0
    assert result_df["Poisson_Enrichment"].notna().all()


def test_del_denoise_with_disynthon_non_unified(disk_datastore):
    """use_disynthon_pairs=True with non_unified strategy produces disynthons + z-score enrichment."""
    config.set_datastore(disk_datastore)
    address = _upload_del_csv(disk_datastore)

    result_address = del_denoise(
        dataset_address=address,
        output_key="denoised_disynthon_non_unified",
        strategy="non_unified",
        use_disynthon_pairs=True,
        add_hit_labels=True,
        hit_percentile=90.0,
    )

    result_df = disk_datastore.get(result_address)
    assert "disynthons" in result_df.columns
    assert "Target_Enrichment_Score" in result_df.columns
    assert "Control_Enrichment_Score" in result_df.columns
    assert "target_hits" in result_df.columns
    assert "control_hits" in result_df.columns
    assert set(result_df["target_hits"].unique()).issubset({0, 1})


def test_del_denoise_with_disynthon_card_metadata(disk_datastore):
    """DataCard includes disynthon metadata when use_disynthon_pairs=True."""
    config.set_datastore(disk_datastore)
    address = _upload_del_csv(disk_datastore)

    result_address = del_denoise(
        dataset_address=address,
        output_key="denoised_disynthon_card",
        strategy="unified",
        use_disynthon_pairs=True,
    )

    card = disk_datastore.get_card(result_address, kind="data")
    assert card.strategy == "unified"
    assert card.use_disynthon_pairs is True
    assert card.n_disynthons > 0
    assert hasattr(card, "n_failed_smiles")
