"""Unit tests for deepchem_server/utils.py aws backend branches."""
from unittest.mock import MagicMock, patch

import pytest

from deepchem_server.utils import _init_datastore, run_job


def test_init_datastore_aws_returns_s3_datastore(monkeypatch):
    monkeypatch.setenv("AWS_BUCKET", "my-bucket")
    with patch("deepchem_server.utils.S3DataStore") as mock_cls:
        mock_cls.return_value = MagicMock()
        _ = _init_datastore("prof", "proj", backend="aws")

    mock_cls.assert_called_once_with(
        profile_name="prof",
        project_name="proj",
        bucket_name="my-bucket",
    )


def test_init_datastore_aws_raises_without_bucket(monkeypatch):
    monkeypatch.delenv("AWS_BUCKET", raising=False)
    with pytest.raises(ValueError, match="AWS_BUCKET"):
        _init_datastore("prof", "proj", backend="aws")


def test_run_job_aws_submits_to_batch_and_returns_job_id(monkeypatch):
    monkeypatch.setenv("AWS_BUCKET", "b")
    monkeypatch.setenv("AWS_BATCH_CPU_JOB_QUEUE", "q")
    monkeypatch.setenv("AWS_BATCH_GPU_JOB_QUEUE", "gq")
    monkeypatch.setenv("AWS_BATCH_CPU_JOB_DEFINITION", "d")
    monkeypatch.setenv("AWS_BATCH_GPU_JOB_DEFINITION", "gd")

    with patch("deepchem_server.aws.batch.boto3") as mock_boto3:
        mock_client = MagicMock()
        mock_boto3.client.return_value = mock_client
        mock_client.submit_job.return_value = {"jobId": "batch-job-123"}

        result = run_job(
            profile_name="prof",
            project_name="proj",
            program={"program_name": "featurize"},
            backend="aws",
        )

    assert result == {"job_id": "batch-job-123", "status": "SUBMITTED"}


def test_run_job_aws_does_not_call_compute_workflow(monkeypatch):
    monkeypatch.setenv("AWS_BUCKET", "b")
    monkeypatch.setenv("AWS_BATCH_CPU_JOB_QUEUE", "q")
    monkeypatch.setenv("AWS_BATCH_GPU_JOB_QUEUE", "gq")
    monkeypatch.setenv("AWS_BATCH_CPU_JOB_DEFINITION", "d")
    monkeypatch.setenv("AWS_BATCH_GPU_JOB_DEFINITION", "gd")

    with patch("deepchem_server.aws.batch.boto3") as mock_boto3, \
         patch("deepchem_server.utils.ComputeWorkflow") as mock_wf:
        mock_client = MagicMock()
        mock_boto3.client.return_value = mock_client
        mock_client.submit_job.return_value = {"jobId": "x"}

        run_job("prof", "proj", {"program_name": "featurize"}, backend="aws")

    mock_wf.assert_not_called()
