"""
Pytest configuration and shared fixtures for pyds tests.
"""

import os
from pathlib import Path
import tempfile

import pytest
import responses

from pyds import BaseClient, Data, Primitives, Settings

@pytest.fixture
def temp_settings_file():
    """Create a temporary settings file for testing."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        temp_file = f.name
    yield temp_file
    # Cleanup
    if os.path.exists(temp_file):
        os.unlink(temp_file)

@pytest.fixture
def settings_instance(temp_settings_file):
    """Create a Settings instance with temporary file."""
    return Settings(
        settings_file=temp_settings_file,
        profile="test-profile",
        project="test-project",
        base_url="http://localhost:8000",
    )

@pytest.fixture
def base_client_instance(settings_instance):
    """Create a BaseClient instance with test settings."""
    return BaseClient(settings=settings_instance)

@pytest.fixture
def data_client_instance(settings_instance):
    """Create a Data client instance with test settings."""
    return Data(settings=settings_instance)

@pytest.fixture
def primitives_client_instance(settings_instance):
    """Create a Primitives client instance with test settings."""
    return Primitives(settings=settings_instance)

@pytest.fixture
def mock_responses():
    """Mock HTTP responses for API calls."""
    with responses.RequestsMock() as rsps:
        yield rsps

@pytest.fixture
def sample_csv_file():
    """Create a sample CSV file for testing uploads."""
    content = """smiles,logp
CC(C)C,1.5
CCO,0.2
CCC,1.0
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(content)
        temp_file = f.name
    yield temp_file
    # Cleanup
    if os.path.exists(temp_file):
        os.unlink(temp_file)

@pytest.fixture
def sample_sdf_file():
    """Create a sample SDF file for testing uploads."""
    content = """
  Mrv2014 01012021

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".sdf", delete=False) as f:
        f.write(content)
        temp_file = f.name
    yield temp_file
    # Cleanup
    if os.path.exists(temp_file):
        os.unlink(temp_file)

@pytest.fixture
def mock_server_response():
    """Standard mock server response for successful operations."""
    return {
        "status": "success",
        "message": "Operation completed successfully",
        "dataset_address": "deepchem://test-profile/test-project/test-file.csv",
    }

@pytest.fixture
def mock_error_response():
    """Standard mock server response for error cases."""
    return {"detail": "Test error message"}

@pytest.fixture
def test_assets_dir():
    """Path to test assets directory."""
    return Path(__file__).parent / "assets"

@pytest.fixture(autouse=True)
def cleanup_settings_files():
    """Automatically cleanup any .pyds.settings.json files created during tests."""
    yield
    # Cleanup any .pyds.settings.json files in the current directory
    for file in [".pyds.settings.json", "test_settings.json"]:
        if os.path.exists(file):
            os.unlink(file)
