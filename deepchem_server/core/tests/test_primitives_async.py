"""Verify primitive routes return AWS Batch job_id dict when run_job returns one."""
from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient

from deepchem_server.main import app

client = TestClient(app)

_BATCH_RESULT = {"job_id": "batch-abc", "status": "SUBMITTED"}


def test_featurize_forwards_batch_result():
    with patch("deepchem_server.routers.primitives.run_job", return_value=_BATCH_RESULT):
        resp = client.post("/primitive/featurize", json={
            "profile_name": "p", "project_name": "q",
            "dataset_address": "deepchem://p/q/data.csv",
            "featurizer": "ecfp", "output": "out", "dataset_column": "smiles",
        })
    assert resp.status_code == 200
    assert resp.json() == _BATCH_RESULT


def test_train_forwards_batch_result():
    with patch("deepchem_server.routers.primitives.run_job", return_value=_BATCH_RESULT):
        resp = client.post("/primitive/train", json={
            "profile_name": "p", "project_name": "q",
            "dataset_address": "deepchem://p/q/data.csv",
            "model_type": "linear_regression", "model_name": "my_model",
        })
    assert resp.status_code == 200
    assert resp.json() == _BATCH_RESULT


def test_del_denoise_forwards_batch_result():
    with patch("deepchem_server.routers.primitives.run_job", return_value=_BATCH_RESULT):
        resp = client.post("/primitive/del/denoise", json={
            "profile_name": "p", "project_name": "q",
            "dataset_address": "deepchem://p/q/del.csv",
            "output_key": "out",
        })
    assert resp.status_code == 200
    assert resp.json() == _BATCH_RESULT
