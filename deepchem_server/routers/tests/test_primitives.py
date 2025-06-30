import pytest
import json
from unittest.mock import Mock, patch
from fastapi.testclient import TestClient
from deepchem_server.main import app


@pytest.fixture
def client():
    """Create a test client for the FastAPI app"""
    return TestClient(app)


class TestFeaturize:
    """Test cases for POST /primitive/featurize endpoint"""

    @patch("deepchem_server.routers.primitives.run_job")
    @patch("deepchem_server.routers.primitives.featurizer_map")
    def test_featurize_success(self, mock_featurizer_map, mock_run_job, client):
        """Test successful featurization request"""
        mock_featurizer_map.keys.return_value = ["ecfp", "graphconv", "weave", "molgraphconv"]

        mock_run_job.return_value = "deepchem://testprofile/testproject/featurized_output.pkl"

        response = client.post(
            "/primitive/featurize",
            json={
                "profile_name": "testprofile",
                "project_name": "testproject",
                "dataset_address": "deepchem://testprofile/testproject/input.csv",
                "featurizer": "ecfp",
                "output": "featurized_output",
                "dataset_column": "smiles",
                "feat_kwargs": {"feat_kwargs": {"radius": 2, "size": 2048}},
                "label_column": "target",
            },
        )

        assert response.status_code == 200
        data = response.json()
        assert data == {
            "featurized_file_address": "deepchem://testprofile/testproject/featurized_output.pkl"
        }

        mock_run_job.assert_called_once()
        call_args = mock_run_job.call_args

        assert call_args[1]["profile_name"] == "testprofile"
        assert call_args[1]["project_name"] == "testproject"

        program = call_args[1]["program"]
        assert program["program_name"] == "featurize"
        assert program["dataset_address"] == "deepchem://testprofile/testproject/input.csv"
        assert program["featurizer"] == "ecfp"
        assert program["output"] == "featurized_output"
        assert program["dataset_column"] == "smiles"
        assert program["label_column"] == "target"
        assert program["feat_kwargs"] == {"radius": 2, "size": 2048}

    @patch("deepchem_server.routers.primitives.featurizer_map")
    def test_featurize_invalid_featurizer(self, mock_featurizer_map, client):
        """Test featurization with invalid featurizer"""
        mock_featurizer_map.keys.return_value = ["ecfp", "graphconv", "weave", "molgraphconv"]

        response = client.post(
            "/primitive/featurize",
            json={
                "profile_name": "testprofile",
                "project_name": "testproject",
                "dataset_address": "deepchem://testprofile/testproject/input.csv",
                "featurizer": "invalid_featurizer",
                "output": "featurized_output",
                "dataset_column": "smiles",
                "feat_kwargs": {"feat_kwargs": {}},
                "label_column": None,
            },
        )

        assert response.status_code == 404
        assert response.json() == {"message": "featurizer invalid"}

    @patch("deepchem_server.routers.primitives.run_job")
    @patch("deepchem_server.routers.primitives.featurizer_map")
    def test_featurize_string_feat_kwargs(self, mock_featurizer_map, mock_run_job, client):
        """Test featurization with string feat_kwargs that need JSON parsing"""
        mock_featurizer_map.keys.return_value = ["ecfp"]
        mock_run_job.return_value = "deepchem://test/test/output.pkl"

        feat_kwargs_str = '{"radius": 3, "size": 1024}'

        response = client.post(
            "/primitive/featurize",
            json={
                "profile_name": "test",
                "project_name": "test",
                "dataset_address": "deepchem://test/test/input.csv",
                "featurizer": "ecfp",
                "output": "output",
                "dataset_column": "smiles",
                "feat_kwargs": {"feat_kwargs": feat_kwargs_str},
            },
        )

        assert response.status_code == 200

        call_args = mock_run_job.call_args
        program = call_args[1]["program"]
        assert program["feat_kwargs"] == {"radius": 3, "size": 1024}

    @patch("deepchem_server.routers.primitives.run_job")
    @patch("deepchem_server.routers.primitives.featurizer_map")
    def test_featurize_boolean_string_conversion(self, mock_featurizer_map, mock_run_job, client):
        """Test conversion of string boolean values in feat_kwargs"""
        mock_featurizer_map.keys.return_value = ["graphconv"]
        mock_run_job.return_value = "deepchem://test/test/output.pkl"

        response = client.post(
            "/primitive/featurize",
            json={
                "profile_name": "test",
                "project_name": "test",
                "dataset_address": "deepchem://test/test/input.csv",
                "featurizer": "graphconv",
                "output": "output",
                "dataset_column": "smiles",
                "feat_kwargs": {
                    "feat_kwargs": {
                        "use_edges": "true",
                        "use_chirality": "false",
                        "use_partial_charge": "none",
                        "max_atoms": "50",
                    }
                },
            },
        )

        assert response.status_code == 200

        call_args = mock_run_job.call_args
        program = call_args[1]["program"]
        feat_kwargs = program["feat_kwargs"]

        assert feat_kwargs["use_edges"] is True
        assert feat_kwargs["use_chirality"] is False
        assert feat_kwargs["use_partial_charge"] is None
        assert feat_kwargs["max_atoms"] == "50"

    @patch("deepchem_server.routers.primitives.run_job")
    @patch("deepchem_server.routers.primitives.featurizer_map")
    def test_featurize_without_label_column(self, mock_featurizer_map, mock_run_job, client):
        """Test featurization without label column (optional parameter)"""
        mock_featurizer_map.keys.return_value = ["ecfp"]
        mock_run_job.return_value = "deepchem://test/test/output.pkl"

        response = client.post(
            "/primitive/featurize",
            json={
                "profile_name": "test",
                "project_name": "test",
                "dataset_address": "deepchem://test/test/input.csv",
                "featurizer": "ecfp",
                "output": "output",
                "dataset_column": "smiles",
                "feat_kwargs": {"feat_kwargs": {}},
            },
        )

        assert response.status_code == 200

        call_args = mock_run_job.call_args
        program = call_args[1]["program"]
        assert program["label_column"] is None

    @patch("deepchem_server.routers.primitives.run_job")
    @patch("deepchem_server.routers.primitives.featurizer_map")
    def test_featurize_empty_feat_kwargs(self, mock_featurizer_map, mock_run_job, client):
        """Test featurization with empty feat_kwargs"""
        mock_featurizer_map.keys.return_value = ["weave"]
        mock_run_job.return_value = "deepchem://test/test/output.pkl"

        response = client.post(
            "/primitive/featurize",
            json={
                "profile_name": "test",
                "project_name": "test",
                "dataset_address": "deepchem://test/test/input.csv",
                "featurizer": "weave",
                "output": "output",
                "dataset_column": "smiles",
                "feat_kwargs": {"feat_kwargs": {}},
            },
        )

        assert response.status_code == 200

        call_args = mock_run_job.call_args
        program = call_args[1]["program"]
        assert program["feat_kwargs"] == {}

    @patch("deepchem_server.routers.primitives.run_job")
    @patch("deepchem_server.routers.primitives.featurizer_map")
    def test_featurize_default_feat_kwargs(self, mock_featurizer_map, mock_run_job, client):
        """Test featurization with default feat_kwargs parameter"""
        mock_featurizer_map.keys.return_value = ["molgraphconv"]
        mock_run_job.return_value = "deepchem://test/test/output.pkl"

        response = client.post(
            "/primitive/featurize",
            json={
                "profile_name": "test",
                "project_name": "test",
                "dataset_address": "deepchem://test/test/input.csv",
                "featurizer": "molgraphconv",
                "output": "output",
                "dataset_column": "smiles",
            },
        )

        assert response.status_code == 200

    def test_featurize_missing_required_params(self, client):
        """Test featurization with missing required parameters"""
        response = client.post(
            "/primitive/featurize",
            json={
                "profile_name": "test",
                "project_name": "test",
            },
        )

        assert response.status_code == 422

    @patch("deepchem_server.routers.primitives.run_job")
    @patch("deepchem_server.routers.primitives.featurizer_map")
    def test_featurize_all_supported_featurizers(self, mock_featurizer_map, mock_run_job, client):
        """Test that all supported featurizers are accepted"""
        supported_featurizers = ["ecfp", "graphconv", "weave", "molgraphconv"]
        mock_featurizer_map.keys.return_value = supported_featurizers
        mock_run_job.return_value = "deepchem://test/test/output.pkl"

        for featurizer in supported_featurizers:
            response = client.post(
                "/primitive/featurize",
                json={
                    "profile_name": "test",
                    "project_name": "test",
                    "dataset_address": "deepchem://test/test/input.csv",
                    "featurizer": featurizer,
                    "output": f"output_{featurizer}",
                    "dataset_column": "smiles",
                    "feat_kwargs": {"feat_kwargs": {}},
                },
            )

            assert response.status_code == 200, f"Failed for featurizer: {featurizer}"

    @patch("deepchem_server.routers.primitives.run_job")
    @patch("deepchem_server.routers.primitives.featurizer_map")
    def test_featurize_complex_feat_kwargs_parsing(self, mock_featurizer_map, mock_run_job, client):
        """Test parsing of complex feat_kwargs structures"""
        mock_featurizer_map.keys.return_value = ["ecfp"]
        mock_run_job.return_value = "deepchem://test/test/output.pkl"

        complex_kwargs = {
            "radius": 2,
            "size": 2048,
            "features": ["true", "false", "none"],
            "nested": {"key": "value"},
            "number": 42,
        }

        response = client.post(
            "/primitive/featurize",
            json={
                "profile_name": "test",
                "project_name": "test",
                "dataset_address": "deepchem://test/test/input.csv",
                "featurizer": "ecfp",
                "output": "output",
                "dataset_column": "smiles",
                "feat_kwargs": {"feat_kwargs": complex_kwargs},
            },
        )

        assert response.status_code == 200

        call_args = mock_run_job.call_args
        program = call_args[1]["program"]
        feat_kwargs = program["feat_kwargs"]

        assert feat_kwargs["radius"] == 2
        assert feat_kwargs["size"] == 2048
        assert feat_kwargs["features"] == [True, False, None]
        assert feat_kwargs["nested"] == {"key": "value"}
        assert feat_kwargs["number"] == 42
