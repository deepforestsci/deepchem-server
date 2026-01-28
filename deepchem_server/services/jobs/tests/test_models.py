"""Tests for Job models."""
import json
import os
from datetime import datetime, timedelta

import pytest

# Set test database before importing models
os.environ["JOBS_DATABASE_URL"] = "sqlite:///:memory:"

from deepchem_server.services.jobs.models import (  # noqa: E402
    Job, JobStatus, get_session, init_db,
)


@pytest.fixture(autouse=True)
def setup_db():
    """Initialize test database."""
    init_db("sqlite:///:memory:")


class TestJobModel:
    """Tests for Job model."""

    def test_create_job(self):
        """Test creating a job."""
        session = get_session()
        try:
            job = Job(
                id="test-job-1",
                program_name="featurize",
                program_json=json.dumps({"program_name": "featurize"}),
                profile="test_profile",
                project="test_project",
            )
            session.add(job)
            session.commit()

            assert job.id == "test-job-1"
            assert job.status == JobStatus.PENDING
            assert job.retry_count == 0
            assert job.max_retries == 3
        finally:
            session.close()

    def test_job_status_transitions(self):
        """Test status transitions."""
        session = get_session()
        try:
            job = Job(
                id="test-job-2",
                program_name="train",
                program_json="{}",
                profile="p",
                project="p",
            )
            session.add(job)
            session.commit()

            # Transition to QUEUED
            job.status = JobStatus.QUEUED
            session.commit()
            assert job.status == JobStatus.QUEUED

            # Transition to RUNNING
            job.status = JobStatus.RUNNING
            job.started_at = datetime.utcnow()
            session.commit()
            assert job.status == JobStatus.RUNNING
            assert job.started_at is not None

            # Transition to COMPLETED
            job.status = JobStatus.COMPLETED
            job.completed_at = datetime.utcnow()
            job.result = json.dumps({"address": "deepchem://p/p/result"})
            session.commit()
            assert job.status == JobStatus.COMPLETED
        finally:
            session.close()

    def test_job_to_dict(self):
        """Test serialization to dict."""
        session = get_session()
        try:
            job = Job(
                id="test-job-3",
                program_name="infer",
                program_json="{}",
                profile="profile",
                project="project",
            )
            session.add(job)
            session.commit()

            d = job.to_dict()
            assert d["id"] == "test-job-3"
            assert d["status"] == "pending"
            assert d["program_name"] == "infer"
            assert "created_at" in d
        finally:
            session.close()

    def test_job_expiry(self):
        """Test expiry date is set correctly."""
        job = Job(
            id="test-job-4",
            program_name="test",
            program_json="{}",
            profile="p",
            project="p",
            retention_days=7,
        )

        # Expiry should be ~7 days from now
        expected = datetime.utcnow() + timedelta(days=7)
        assert abs((job.expires_at - expected).total_seconds()) < 5
