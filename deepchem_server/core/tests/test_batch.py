"""Unit tests for deepchem_server/aws/batch.py."""
import json
from unittest.mock import MagicMock, patch

import pytest


@pytest.fixture
def batch_env(monkeypatch):
    monkeypatch.setenv("AWS_BUCKET", "test-bucket")
    monkeypatch.setenv("AWS_BATCH_CPU_JOB_QUEUE", "cpu-queue")
    monkeypatch.setenv("AWS_BATCH_GPU_JOB_QUEUE", "gpu-queue")
    monkeypatch.setenv("AWS_BATCH_CPU_JOB_DEFINITION", "cpu-job-def")
    monkeypatch.setenv("AWS_BATCH_GPU_JOB_DEFINITION", "gpu-job-def")


def _mock_boto3_batch(job_id="job-abc"):
    mock_client = MagicMock()
    mock_client.submit_job.return_value = {"jobId": job_id}
    mock_boto3 = MagicMock()
    mock_boto3.client.return_value = mock_client
    return mock_boto3, mock_client


def test_submit_cpu_primitive_uses_cpu_queue(batch_env):
    mock_boto3, mock_client = _mock_boto3_batch()
    with patch("deepchem_server.aws.batch.boto3", mock_boto3):
        from deepchem_server.aws.batch import submit_job
        job_id = submit_job({"program_name": "featurize"}, "prof", "proj")

    assert job_id == "job-abc"
    kwargs = mock_client.submit_job.call_args[1]
    assert kwargs["jobQueue"] == "cpu-queue"
    assert kwargs["jobDefinition"] == "cpu-job-def"
    reqs = {r["type"]: r["value"] for r in kwargs["containerOverrides"]["resourceRequirements"]}
    assert "GPU" not in reqs


def test_submit_gpu_primitive_uses_gpu_queue(batch_env):
    mock_boto3, mock_client = _mock_boto3_batch("job-gpu")
    with patch("deepchem_server.aws.batch.boto3", mock_boto3):
        from deepchem_server.aws.batch import submit_job
        job_id = submit_job({"program_name": "train"}, "prof", "proj")

    assert job_id == "job-gpu"
    kwargs = mock_client.submit_job.call_args[1]
    assert kwargs["jobQueue"] == "gpu-queue"
    assert kwargs["jobDefinition"] == "gpu-job-def"
    reqs = {r["type"]: r["value"] for r in kwargs["containerOverrides"]["resourceRequirements"]}
    assert reqs["GPU"] == "1"


def test_submit_unknown_primitive_falls_back_to_cpu_defaults(batch_env):
    mock_boto3, mock_client = _mock_boto3_batch()
    with patch("deepchem_server.aws.batch.boto3", mock_boto3):
        from deepchem_server.aws.batch import submit_job
        submit_job({"program_name": "unknown_op"}, "prof", "proj")

    kwargs = mock_client.submit_job.call_args[1]
    assert kwargs["jobQueue"] == "cpu-queue"
    reqs = {r["type"]: r["value"] for r in kwargs["containerOverrides"]["resourceRequirements"]}
    assert reqs["VCPU"] == "4"
    assert reqs["MEMORY"] == "8192"


def test_submit_embeds_program_profile_project_bucket_in_command(batch_env):
    program = {"program_name": "del_denoise", "dataset_address": "deepchem://a/b/c.csv"}
    mock_boto3, mock_client = _mock_boto3_batch()
    with patch("deepchem_server.aws.batch.boto3", mock_boto3):
        from deepchem_server.aws.batch import submit_job
        submit_job(program, "my-profile", "my-project")

    cmd = mock_client.submit_job.call_args[1]["containerOverrides"]["command"]
    assert "--program" in cmd
    assert "--profile" in cmd
    assert "--project" in cmd
    assert "--bucket" in cmd
    assert json.loads(cmd[cmd.index("--program") + 1]) == program
    assert cmd[cmd.index("--profile") + 1] == "my-profile"
    assert cmd[cmd.index("--project") + 1] == "my-project"
    assert cmd[cmd.index("--bucket") + 1] == "test-bucket"


def test_get_job_status_succeeded():
    mock_boto3 = MagicMock()
    mock_client = MagicMock()
    mock_boto3.client.return_value = mock_client
    mock_client.describe_jobs.return_value = {
        "jobs": [{
            "status": "SUCCEEDED",
            "statusReason": "Essential container in task exited",
            "container": {"logStreamName": "group/stream/123"},
        }]
    }
    with patch("deepchem_server.aws.batch.boto3", mock_boto3):
        from deepchem_server.aws.batch import get_job_status
        result = get_job_status("job-abc")

    assert result["job_id"] == "job-abc"
    assert result["status"] == "SUCCEEDED"
    assert result["log_stream"] == "group/stream/123"


def test_get_job_status_not_found():
    mock_boto3 = MagicMock()
    mock_client = MagicMock()
    mock_boto3.client.return_value = mock_client
    mock_client.describe_jobs.return_value = {"jobs": []}
    with patch("deepchem_server.aws.batch.boto3", mock_boto3):
        from deepchem_server.aws.batch import get_job_status
        result = get_job_status("ghost-id")

    assert result["status"] == "NOT_FOUND"
