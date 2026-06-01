import json
import os
import uuid

import boto3

RESOURCE_PROFILES = {
    "featurize":                    {"vcpu": "4",  "memory": "8192",  "gpu": False},
    "train":                        {"vcpu": "8",  "memory": "32768", "gpu": True},
    "infer":                        {"vcpu": "4",  "memory": "16384", "gpu": True},
    "evaluate":                     {"vcpu": "4",  "memory": "16384", "gpu": True},
    "train_valid_test_split":       {"vcpu": "4",  "memory": "8192",  "gpu": False},
    "partition":                    {"vcpu": "4",  "memory": "8192",  "gpu": False},
    "generate_pose":                {"vcpu": "4",  "memory": "8192",  "gpu": False},
    "relative_binding_free_energy": {"vcpu": "16", "memory": "65536", "gpu": False},
    "collate_rbfe_results":         {"vcpu": "4",  "memory": "8192",  "gpu": False},
    "del_denoise":                  {"vcpu": "4",  "memory": "8192",  "gpu": False},
}
_DEFAULT = {"vcpu": "4", "memory": "8192", "gpu": False}


def submit_job(program: dict, profile_name: str, project_name: str) -> str:
    """Submit a primitive program to AWS Batch and return the Batch job ID."""
    profile  = RESOURCE_PROFILES.get(program.get("program_name", ""), _DEFAULT)
    use_gpu  = profile["gpu"]

    job_queue   = os.environ["AWS_BATCH_GPU_JOB_QUEUE"       if use_gpu else "AWS_BATCH_CPU_JOB_QUEUE"]
    job_def     = os.environ["AWS_BATCH_GPU_JOB_DEFINITION"  if use_gpu else "AWS_BATCH_CPU_JOB_DEFINITION"]
    bucket      = os.environ["AWS_BUCKET"]

    resource_reqs = [
        {"type": "VCPU",   "value": profile["vcpu"]},
        {"type": "MEMORY", "value": profile["memory"]},
    ]
    if use_gpu:
        resource_reqs.append({"type": "GPU", "value": "1"})

    job_name = f"{program.get('program_name', 'job')}-{uuid.uuid4().hex[:8]}"

    response = boto3.client("batch").submit_job(
        jobName=job_name,
        jobQueue=job_queue,
        jobDefinition=job_def,
        containerOverrides={
            "command": [
                "python", "-m", "deepchem_server.aws.command",
                "--program", json.dumps(program),
                "--profile", profile_name,
                "--project", project_name,
                "--bucket",  bucket,
            ],
            "resourceRequirements": resource_reqs,
        },
    )
    return response["jobId"]


def get_job_status(job_id: str) -> dict:
    """Return status dict for a Batch job ID."""
    jobs = boto3.client("batch").describe_jobs(jobs=[job_id])["jobs"]
    if not jobs:
        return {"job_id": job_id, "status": "NOT_FOUND", "status_reason": "", "log_stream": ""}
    job = jobs[0]
    return {
        "job_id":        job_id,
        "status":        job["status"],
        "status_reason": job.get("statusReason", ""),
        "log_stream":    job.get("container", {}).get("logStreamName", ""),
    }
