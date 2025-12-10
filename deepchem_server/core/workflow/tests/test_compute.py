import pandas as pd
import pytest

from deepchem_server.core.common import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.workflow.compute import ComputeWorkflow, program_map


def test_program_map_contains_expected_programs():
    """Test that program_map contains all expected programs."""
    expected_programs = [
        'featurize',
        'train',
        'evaluate',
        'infer',
        'train_valid_test_split',
        'generate_pose',
        'relative_binding_free_energy',
        'collate_rbfe_results',
    ]
    for program in expected_programs:
        assert program in program_map, f"Expected program '{program}' not in program_map"


def test_compute_workflow_featurize(disk_datastore):
    """Test ComputeWorkflow with featurize program."""
    config.set_datastore(disk_datastore)

    # Create test data
    df = pd.DataFrame(["CCC", "CCCCC"], columns=["smiles"])
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    data_address = disk_datastore.upload_data_from_memory(df, "test.csv", card)

    # Create workflow program
    program = {
        'program_name': 'featurize',
        'dataset_address': data_address,
        'featurizer': 'ecfp',
        'output': 'feat_test',
        'dataset_column': 'smiles',
    }

    workflow = ComputeWorkflow(program)
    result = workflow.execute()

    assert result is not None
    assert result.startswith('deepchem://')

    # Verify featurized data
    data = disk_datastore.get(result)
    assert data.X.shape == (2, 2048)


def test_compute_workflow_missing_program_name():
    """Test that ComputeWorkflow raises error when program_name is missing."""
    program = {
        'dataset_address': 'deepchem://test/user/data.csv',
        'featurizer': 'ecfp',
    }

    workflow = ComputeWorkflow(program)

    with pytest.raises(ValueError, match="program_name not found in program"):
        workflow.execute()


def test_compute_workflow_invalid_program_name():
    """Test that ComputeWorkflow raises error for invalid program_name."""
    program = {
        'program_name': 'invalid_program',
        'dataset_address': 'deepchem://test/user/data.csv',
    }

    workflow = ComputeWorkflow(program)

    with pytest.raises(ValueError, match="invalid_program not in available programs"):
        workflow.execute()


def test_compute_workflow_train(disk_datastore):
    """Test ComputeWorkflow with train program."""
    config.set_datastore(disk_datastore)

    # Create and featurize test data
    df = pd.DataFrame([["CCC", 0], ["CCCCC", 1]], columns=["smiles", "label"])
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    data_address = disk_datastore.upload_data_from_memory(df, "test.csv", card)

    # Featurize first
    feat_program = {
        'program_name': 'featurize',
        'dataset_address': data_address,
        'featurizer': 'ecfp',
        'output': 'feat_test',
        'dataset_column': 'smiles',
        'label_column': 'label',
    }
    feat_workflow = ComputeWorkflow(feat_program)
    feat_address = feat_workflow.execute()

    # Train model
    train_program = {
        'program_name': 'train',
        'model_type': 'random_forest_classifier',
        'dataset_address': feat_address,
        'model_name': 'test_model',
        'init_kwargs': {},
        'train_kwargs': {},
    }

    train_workflow = ComputeWorkflow(train_program)
    model_address = train_workflow.execute()

    assert model_address is not None
    assert model_address.startswith('deepchem://')

    # Verify model was saved
    model = disk_datastore.get_model(model_address)
    assert model is not None
