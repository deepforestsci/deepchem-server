import pandas as pd
import pytest

from deepchem_server.core import config, merge
from deepchem_server.core.cards import DataCard


def test_merge_two_csv_files(disk_datastore):
    """Test merging two CSV files."""
    config.set_datastore(disk_datastore)

    df1 = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
    df2 = pd.DataFrame({"col1": [4, 5], "col2": ["d", "e"]})

    card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    address1 = disk_datastore.upload_data_from_memory(df1, "test1.csv", card)
    address2 = disk_datastore.upload_data_from_memory(df2, "test2.csv", card)

    merged_address = merge.merge([address1, address2], "merged_output.csv")
    merged_df = disk_datastore.get(merged_address)

    assert len(merged_df) == 5
    assert list(merged_df.columns) == ["col1", "col2"]
    assert merged_df["col1"].tolist() == [1, 2, 3, 4, 5]
    assert merged_df["col2"].tolist() == ["a", "b", "c", "d", "e"]
    assert merged_address == "deepchem://test/user/merged_output.csv"


def test_merge_multiple_csv_files(disk_datastore):
    """Test merging more than two CSV files."""
    config.set_datastore(disk_datastore)

    df1 = pd.DataFrame({"x": [1, 2], "y": [10, 20]})
    df2 = pd.DataFrame({"x": [3], "y": [30]})
    df3 = pd.DataFrame({"x": [4, 5, 6], "y": [40, 50, 60]})

    card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    address1 = disk_datastore.upload_data_from_memory(df1, "file1.csv", card)
    address2 = disk_datastore.upload_data_from_memory(df2, "file2.csv", card)
    address3 = disk_datastore.upload_data_from_memory(df3, "file3.csv", card)

    merged_address = merge.merge([address1, address2, address3], "merged_multi.csv")

    merged_df = disk_datastore.get(merged_address)
    assert len(merged_df) == 6
    assert merged_df["x"].tolist() == [1, 2, 3, 4, 5, 6]
    assert merged_df["y"].tolist() == [10, 20, 30, 40, 50, 60]


def test_merge_with_string_input(disk_datastore):
    """Test merge with string input for dataset_addresses (tests ast.literal_eval path)."""
    config.set_datastore(disk_datastore)

    df1 = pd.DataFrame({"name": ["Alice", "Bob"], "age": [25, 30]})
    df2 = pd.DataFrame({"name": ["Charlie"], "age": [35]})

    card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    address1 = disk_datastore.upload_data_from_memory(df1, "people1.csv", card)
    address2 = disk_datastore.upload_data_from_memory(df2, "people2.csv", card)

    addresses_str = f"['{address1}', '{address2}']"
    merged_address = merge.merge(addresses_str, "merged_people.csv")

    merged_df = disk_datastore.get(merged_address)
    assert len(merged_df) == 3
    assert merged_df["name"].tolist() == ["Alice", "Bob", "Charlie"]
    assert merged_df["age"].tolist() == [25, 30, 35]


def test_merge_nested_path(disk_datastore):
    """Test merging with nested folder paths."""
    config.set_datastore(disk_datastore)

    df1 = pd.DataFrame({"value": [100, 200]})
    df2 = pd.DataFrame({"value": [300]})

    card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    address1 = disk_datastore.upload_data_from_memory(df1, "input/data1.csv", card)
    address2 = disk_datastore.upload_data_from_memory(df2, "input/data2.csv", card)

    merged_address = merge.merge([address1, address2], "output/merged_data.csv")

    merged_df = disk_datastore.get(merged_address)
    assert len(merged_df) == 3
    assert merged_df["value"].tolist() == [100, 200, 300]
    assert merged_address == "deepchem://test/user/output/merged_data.csv"


def test_merge_preserves_column_order(disk_datastore):
    """Test that merge preserves the column order from the first file."""
    config.set_datastore(disk_datastore)

    df1 = pd.DataFrame({"z": [1], "a": [2], "m": [3]})
    df2 = pd.DataFrame({"z": [4], "a": [5], "m": [6]})

    card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    address1 = disk_datastore.upload_data_from_memory(df1, "cols1.csv", card)
    address2 = disk_datastore.upload_data_from_memory(df2, "cols2.csv", card)

    merged_address = merge.merge([address1, address2], "merged_cols.csv")

    merged_df = disk_datastore.get(merged_address)
    assert list(merged_df.columns) == ["z", "a", "m"]
    assert len(merged_df) == 2


def test_merge_empty_dataframes(disk_datastore):
    """Test merging CSV files where some are empty."""
    config.set_datastore(disk_datastore)

    df1 = pd.DataFrame({"col": []})
    df2 = pd.DataFrame({"col": [1, 2, 3]})

    card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    address1 = disk_datastore.upload_data_from_memory(df1, "empty.csv", card)
    address2 = disk_datastore.upload_data_from_memory(df2, "nonempty.csv", card)

    merged_address = merge.merge([address1, address2], "merged_empty.csv")

    merged_df = disk_datastore.get(merged_address)
    assert len(merged_df) == 3
    assert merged_df["col"].tolist() == [1, 2, 3]


def test_merge_non_csv_raises_error(disk_datastore):
    """Test that merging non-CSV files raises a TypeError."""
    config.set_datastore(disk_datastore)

    df1 = pd.DataFrame({"col": [1, 2]})
    df2 = pd.DataFrame({"col": [3, 4]})

    card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    address1 = disk_datastore.upload_data_from_memory(df1, "test1.csv", card)
    address2 = disk_datastore.upload_data_from_memory(df2, "test2.csv", card)

    with pytest.raises(TypeError, match="merge operation is only available for csv file types"):
        merge.merge([address1, address2], "output.json")


def test_merge_single_file(disk_datastore):
    """Test merging a single CSV file (edge case)."""
    config.set_datastore(disk_datastore)

    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})

    card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    address = disk_datastore.upload_data_from_memory(df, "single.csv", card)

    merged_address = merge.merge([address], "merged_single.csv")
    merged_df = disk_datastore.get(merged_address)

    assert len(merged_df) == 3
    assert merged_df["a"].tolist() == [1, 2, 3]
    assert merged_df["b"].tolist() == [4, 5, 6]


def test_merge_with_missing_datastore(disk_datastore):
    """Test that merge raises ValueError when datastore is not set."""
    config.set_datastore(None)

    df1 = pd.DataFrame({"col": [1, 2]})
    df2 = pd.DataFrame({"col": [3, 4]})

    card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    address1 = disk_datastore.upload_data_from_memory(df1, "test1.csv", card)
    address2 = disk_datastore.upload_data_from_memory(df2, "test2.csv", card)

    with pytest.raises(ValueError, match="Datastore not set"):
        merge.merge([address1, address2], "output.csv")


def test_merge_with_different_column_types(disk_datastore):
    """Test merging CSV files with different data types."""
    config.set_datastore(disk_datastore)

    df1 = pd.DataFrame({"id": [1, 2], "name": ["a", "b"], "value": [1.5, 2.5]})
    df2 = pd.DataFrame({"id": [3], "name": ["c"], "value": [3.5]})

    card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    address1 = disk_datastore.upload_data_from_memory(df1, "types1.csv", card)
    address2 = disk_datastore.upload_data_from_memory(df2, "types2.csv", card)

    merged_address = merge.merge([address1, address2], "merged_types.csv")

    merged_df = disk_datastore.get(merged_address)
    assert len(merged_df) == 3
    assert merged_df["name"].tolist() == ["a", "b", "c"]
    assert len(merged_df["value"]) == 3
