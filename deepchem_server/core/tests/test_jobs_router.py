"""Unit tests for GET /jobs/{job_id}."""
from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient

from deepchem_server.main import app

client = TestClient(app)


def test_get_job_status_returns_succeeded():
    with patch("deepchem_server.routers.jobs.get_job_status") as mock_fn:
        mock_fn.return_value = {
            "job_id":        "abc-123",
            "status":        "SUCCEEDED",
            "status_reason": "",
            "log_stream":    "group/stream/1",
        }
        response = client.get("/jobs/abc-123")

    assert response.status_code == 200
    body = response.json()
    assert body["job_id"] == "abc-123"
    assert body["status"] == "SUCCEEDED"


def test_get_job_status_running():
    with patch("deepchem_server.routers.jobs.get_job_status") as mock_fn:
        mock_fn.return_value = {
            "job_id": "run-1", "status": "RUNNING",
            "status_reason": "", "log_stream": "",
        }
        response = client.get("/jobs/run-1")

    assert response.status_code == 200
    assert response.json()["status"] == "RUNNING"


def test_get_job_status_not_found_returns_404():
    with patch("deepchem_server.routers.jobs.get_job_status") as mock_fn:
        mock_fn.return_value = {
            "job_id": "ghost", "status": "NOT_FOUND",
            "status_reason": "", "log_stream": "",
        }
        response = client.get("/jobs/ghost")

    assert response.status_code == 404
