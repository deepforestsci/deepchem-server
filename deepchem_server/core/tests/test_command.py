"""Unit tests for deepchem_server/aws/command.py."""
import argparse
import json
from unittest.mock import MagicMock, patch


def _make_args(program_name="featurize", **extra):
    program = {"program_name": program_name, **extra}
    return argparse.Namespace(
        program=json.dumps(program),
        profile="test-profile",
        project="test-project",
        bucket="test-bucket",
    )


def test_command_inits_s3_datastore_with_correct_args():
    args = _make_args()
    with patch("deepchem_server.aws.command.S3DataStore") as mock_ds_cls, \
         patch("deepchem_server.aws.command.config"), \
         patch("deepchem_server.aws.command.ComputeWorkflow") as mock_wf_cls:
        mock_wf_cls.return_value = MagicMock()
        from deepchem_server.aws.command import main
        main(args)

    mock_ds_cls.assert_called_once_with(
        profile_name="test-profile",
        project_name="test-project",
        bucket_name="test-bucket",
    )


def test_command_sets_datastore_on_config():
    args = _make_args()
    with patch("deepchem_server.aws.command.S3DataStore") as mock_ds_cls, \
         patch("deepchem_server.aws.command.config") as mock_config, \
         patch("deepchem_server.aws.command.ComputeWorkflow") as mock_wf_cls:
        mock_ds = MagicMock()
        mock_ds_cls.return_value = mock_ds
        mock_wf_cls.return_value = MagicMock()
        from deepchem_server.aws.command import main
        main(args)

    mock_config.set_datastore.assert_called_once_with(mock_ds)


def test_command_runs_correct_workflow():
    program = {"program_name": "del_denoise", "dataset_address": "deepchem://a/b/c.csv"}
    args = argparse.Namespace(
        program=json.dumps(program),
        profile="p",
        project="q",
        bucket="b",
    )
    with patch("deepchem_server.aws.command.S3DataStore"), \
         patch("deepchem_server.aws.command.config"), \
         patch("deepchem_server.aws.command.ComputeWorkflow") as mock_wf_cls:
        mock_wf = MagicMock()
        mock_wf_cls.return_value = mock_wf
        from deepchem_server.aws.command import main
        main(args)

    mock_wf_cls.assert_called_once_with(program)
    mock_wf.execute.assert_called_once()
