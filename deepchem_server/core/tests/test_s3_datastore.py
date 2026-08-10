"""Unit tests for S3DataStore using moto to mock AWS S3."""

import os
from unittest.mock import patch

import boto3
import pandas as pd
import pytest
from moto import mock_aws

from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.datastore import S3DataStore


BUCKET = "test-deepchem-bucket"
PROFILE = "testprofile"
PROJECT = "testproject"


@pytest.fixture
def aws_credentials():
    """Mock AWS credentials"""
    os.environ["AWS_ACCESS_KEY_ID"] = "testing"
    os.environ["AWS_SECRET_ACCESS_KEY"] = "testing"
    os.environ["AWS_SECURITY_TOKEN"] = "testing"
    os.environ["AWS_SESSION_TOKEN"] = "testing"
    os.environ["AWS_DEFAULT_REGION"] = "us-east-1"
    yield
    for key in (
            "AWS_ACCESS_KEY_ID",
            "AWS_SECRET_ACCESS_KEY",
            "AWS_SECURITY_TOKEN",
            "AWS_SESSION_TOKEN",
            "AWS_DEFAULT_REGION",
    ):
        os.environ.pop(key, None)


@pytest.fixture
def s3_bucket(aws_credentials):
    """Create a mocked S3 bucket and yield the bucket name."""
    with mock_aws():
        boto3.client("s3", region_name="us-east-1").create_bucket(Bucket=BUCKET)
        yield BUCKET


@pytest.fixture
def s3_datastore(s3_bucket):
    """S3DataStore backed by the mocked bucket."""
    return S3DataStore(profile_name=PROFILE, project_name=PROJECT, bucket_name=s3_bucket)


@pytest.fixture
def tmp_csv(tmp_path):
    path = tmp_path / "data.csv"
    pd.DataFrame({"x": [1, 2, 3], "y": [4, 5, 6]}).to_csv(path, index=False)
    return str(path)


@pytest.fixture
def tmp_parquet(tmp_path):
    path = tmp_path / "data.parquet"
    pd.DataFrame({"a": [10, 20], "b": [30, 40]}).to_parquet(path, index=False)
    return str(path)


@pytest.fixture
def tmp_text(tmp_path):
    path = tmp_path / "note.txt"
    path.write_text("hello deepchem")
    return str(path)


def _csv_card(**kwargs):
    return DataCard(address="", file_type="csv", data_type="pandas.DataFrame", **kwargs)


def _parquet_card(**kwargs):
    return DataCard(address="", file_type="parquet", data_type="pandas.DataFrame", **kwargs)


def _text_card(**kwargs):
    return DataCard(address="", file_type="txt", data_type="text/plain", **kwargs)


def test_s3_datastore_init_fails_without_bucket_name(s3_bucket):
    with pytest.raises(AssertionError):
        S3DataStore(profile_name=PROFILE, project_name=PROJECT, bucket_name="")


def test_upload_csv_returns_deepchem_address(s3_datastore, tmp_csv):
    address = s3_datastore.upload_data("data.csv", tmp_csv, _csv_card())
    assert address == f"deepchem://{PROFILE}/{PROJECT}/data.csv"


def test_upload_csv_and_retrieve_as_dataframe(s3_datastore, tmp_csv):
    address = s3_datastore.upload_data("data.csv", tmp_csv, _csv_card())
    df = s3_datastore.get_data(address)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["x", "y"]
    assert len(df) == 3


def test_upload_csv_card_shape_is_set(s3_datastore, tmp_csv):
    card = _csv_card()
    s3_datastore.upload_data("data.csv", tmp_csv, card)
    assert card.shape == (3, 2)


def test_upload_parquet_and_retrieve_as_dataframe(s3_datastore, tmp_parquet):
    address = s3_datastore.upload_data("data.parquet", tmp_parquet, _parquet_card())
    df = s3_datastore.get_data(address)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["a", "b"]
    assert len(df) == 2


def test_upload_text_file_and_retrieve_as_bytes(s3_datastore, tmp_text):
    address = s3_datastore.upload_data("note.txt", tmp_text, _text_card())
    result = s3_datastore.get_data(address)
    assert isinstance(result, bytes)
    assert b"hello deepchem" in result


def test_get_card_after_upload(s3_datastore, tmp_csv):
    address = s3_datastore.upload_data("data.csv", tmp_csv, _csv_card(description="test csv"))
    card = s3_datastore.get_card(address, kind="data")
    assert isinstance(card, DataCard)
    assert card.file_type == "csv"
    assert card.data_type == "pandas.DataFrame"


def test_get_dispatches_to_data(s3_datastore, tmp_csv):
    address = s3_datastore.upload_data("data.csv", tmp_csv, _csv_card())
    df = s3_datastore.get(address, kind="data")
    assert isinstance(df, pd.DataFrame)


def test_get_with_cdc_suffix_returns_card(s3_datastore, tmp_csv):
    address = s3_datastore.upload_data("data.csv", tmp_csv, _csv_card())
    card = s3_datastore.get(address + ".cdc")
    assert isinstance(card, DataCard)


def test_list_data_shows_uploaded_file(s3_datastore, tmp_csv):
    s3_datastore.upload_data("data.csv", tmp_csv, _csv_card())
    listing = s3_datastore.list_data()
    assert "data.csv" in listing


def test_list_data_empty_for_fresh_datastore(s3_datastore):
    listing = s3_datastore.list_data()
    assert listing == ""


def test_list_data_shows_multiple_files(s3_datastore, tmp_csv, tmp_text):
    s3_datastore.upload_data("alpha.csv", tmp_csv, _csv_card())
    s3_datastore.upload_data("beta.txt", tmp_text, _text_card())
    listing = s3_datastore.list_data()
    assert "alpha.csv" in listing
    assert "beta.txt" in listing


def test_exists_returns_false_for_missing_file(s3_datastore):
    assert not s3_datastore.exists("deepchem://testprofile/testproject/ghost.csv", kind="data")


def test_exists_returns_true_after_upload(s3_datastore, tmp_csv):
    address = s3_datastore.upload_data("data.csv", tmp_csv, _csv_card())
    assert s3_datastore.exists(address, kind="data")


def test_exists_returns_false_after_delete(s3_datastore, tmp_csv):
    address = s3_datastore.upload_data("data.csv", tmp_csv, _csv_card())
    s3_datastore.delete_object(address, kind="data")
    assert not s3_datastore.exists(address, kind="data")


def test_delete_removes_object_and_card_from_s3(s3_datastore, tmp_csv):
    address = s3_datastore.upload_data("data.csv", tmp_csv, _csv_card())
    s3_datastore.delete_object(address, kind="data")
    assert s3_datastore.get_data(address) is None


def test_download_object_writes_file_to_path(s3_datastore, tmp_csv, tmp_path):
    address = s3_datastore.upload_data("data.csv", tmp_csv, _csv_card())
    dest = str(tmp_path / "downloaded.csv")
    s3_datastore.download_object(address, filename=dest)
    assert os.path.isfile(dest)
    df = pd.read_csv(dest)
    assert len(df) == 3


def test_download_object_uses_temp_path_when_no_filename(s3_datastore, tmp_csv):
    address = s3_datastore.upload_data("data.csv", tmp_csv, _csv_card())
    with pytest.warns(FutureWarning):
        path = s3_datastore.download_object(address)
    assert os.path.isfile(str(path))


def test_download_object_raises_for_missing_address(s3_datastore):
    with pytest.raises(ValueError, match="does not exist"):
        s3_datastore.download_object("deepchem://testprofile/testproject/ghost.csv")


def test_get_file_size_matches_original(s3_datastore, tmp_csv):
    address = s3_datastore.upload_data("data.csv", tmp_csv, _csv_card())
    expected = os.path.getsize(tmp_csv)
    assert s3_datastore.get_file_size(address) == expected


def test_add_dir_appears_in_listing(s3_datastore):
    s3_datastore.add_dir("experiments")
    listing = s3_datastore.list_data()
    assert "experiments" in listing


def test_add_dir_raises_if_already_exists(s3_datastore, tmp_csv):
    s3_datastore.upload_data("mydir/data.csv", tmp_csv, _csv_card())
    with pytest.raises(ValueError, match="already exists"):
        s3_datastore.add_dir("mydir")


def test_upload_directory_and_get_dir(s3_datastore, tmp_path):
    d = tmp_path / "subdir"
    d.mkdir()
    (d / "a.txt").write_text("aaa")
    (d / "b.txt").write_text("bbb")
    card = DataCard(address="", file_type="dir", data_type="text/plain")
    address = s3_datastore.upload_data("mydir", str(d), card)

    local_dir = s3_datastore.get_dir(address)
    assert os.path.isdir(local_dir)
    files = os.listdir(local_dir)
    assert "a.txt" in files
    assert "b.txt" in files


def test_get_data_sample_passes_range_header(s3_datastore, tmp_csv):
    address = s3_datastore.upload_data("data.csv", tmp_csv, _csv_card())

    with patch.object(s3_datastore._s3_client, "get_object", wraps=s3_datastore._s3_client.get_object) as spy:
        s3_datastore.get_data(address, fetch_sample=True)
        call_kwargs = spy.call_args.kwargs
        assert "Range" in call_kwargs
        assert call_kwargs["Range"] == s3_datastore.range_limit_header


def test_get_data_no_sample_omits_range_header(s3_datastore, tmp_csv):
    address = s3_datastore.upload_data("data.csv", tmp_csv, _csv_card())

    with patch.object(s3_datastore._s3_client, "get_object", wraps=s3_datastore._s3_client.get_object) as spy:
        s3_datastore.get_data(address, fetch_sample=False)
        call_kwargs = spy.call_args.kwargs
        assert "Range" not in call_kwargs
