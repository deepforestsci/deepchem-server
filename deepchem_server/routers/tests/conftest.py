import pytest
import os
import tempfile
import sys
from unittest.mock import patch, MagicMock

sys.modules["deepchem"] = MagicMock()
sys.modules["deepchem.feat"] = MagicMock()
sys.modules["deepchem.data"] = MagicMock()
sys.modules["rdkit"] = MagicMock()
sys.modules["rdkit.Chem"] = MagicMock()

from fastapi.testclient import TestClient
from deepchem_server.main import app


@pytest.fixture(scope="session")
def test_client():
    """Create a test client for the FastAPI app that can be shared across test sessions"""
    return TestClient(app)


@pytest.fixture(autouse=True)
def mock_environment():
    """Mock environment variables for testing"""
    with patch.dict(os.environ, {"DATADIR": "/tmp/test_datastore", "ENVIRONMENT": "test"}):
        yield


@pytest.fixture
def temp_datastore():
    """Create a temporary datastore directory for tests"""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def sample_csv_content():
    """Sample CSV content for testing uploads"""
    return b"smiles,label\nCCC,1\nCCCC,0\nCCCCC,1"


FILE_TYPE_TEST_DATA = [
    ("test.csv", "pandas.DataFrame"),
    ("test.parquet", "pandas.DataFrame"),
    ("test.pdb", "text/plain"),
    ("test.sdf", "text/plain"),
    ("test.fasta", "text/plain"),
    ("test.fastq", "text/plain"),
    ("test.txt", "text/plain"),
    ("test.xml", "text/plain"),
    ("test.pdbqt", "text/plain"),
    ("test.smi", "text/plain"),
    ("test.smiles", "text/plain"),
    ("test.cxsmiles", "text/plain"),
    ("test.json", "text/plain"),
    ("test.dcd", "binary"),
    ("test.bz2", "binary"),
    ("test.zip", "binary"),
    ("test.onnx", "binary"),
    ("test.hdf5", "binary"),
    ("test.png", "png"),
    ("test.unknown", ""),
]


@pytest.fixture(params=FILE_TYPE_TEST_DATA)
def file_type_data(request):
    """Parametrized fixture for testing different file types"""
    return request.param


SUPPORTED_FEATURIZERS = ["ecfp", "graphconv", "weave", "molgraphconv"]


@pytest.fixture(params=SUPPORTED_FEATURIZERS)
def featurizer_name(request):
    """Parametrized fixture for testing different featurizers"""
    return request.param


@pytest.fixture
def mock_successful_upload():
    """Mock successful upload response"""
    return "deepchem://testprofile/testproject/test_file.csv"


@pytest.fixture
def mock_datastore_objects():
    """Mock datastore objects for listing and searching"""
    return [
        "dataset1.csv",
        "dataset2.sdf",
        "experiment_data.txt",
        "molecular_features.pkl",
        "model_weights.h5",
    ]
