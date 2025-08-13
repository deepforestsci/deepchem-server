# PyDS Tests

This directory contains comprehensive tests for the pyds package using pytest.

## Test Structure

The tests are organized into separate files for each module:

- `conftest.py` - Shared fixtures and test configuration
- `test_settings.py` - Tests for the Settings module
- `test_base_client.py` - Tests for the BaseClient module  
- `test_data.py` - Tests for the Data module
- `test_primitives.py` - Tests for the Primitives module (organized by primitive operation)
- `test_init.py` - Tests for package initialization and imports

### Primitives Test Organization

The `test_primitives.py` file is organized into focused test classes:

- `TestPrimitivesCommonUtils` - Common utilities, initialization, and error handling tests
- `TestFeaturize` - Tests for featurization operations
- `TestTrain` - Tests for training operations  
- `TestEvaluate` - Tests for evaluation operations
- `TestInfer` - Tests for inference operations

## Running Tests

### Install Dependencies

First, install the testing dependencies:

```bash
pip install -r requirements.txt
```

### Run All Tests

To run all tests:

```bash
pytest
```

### Run Unit Tests Only

To run only unit tests (excluding integration tests):

```bash
pytest -m "not integration"
```

### Run Integration Tests

To run integration tests (requires a running DeepChem server):

```bash
pytest -m integration
```

### Run Tests with Coverage

To run tests with coverage report:

```bash
pytest --cov=pyds --cov-report=html
```

### Run Specific Test Files

To run tests from a specific file:

```bash
pytest tests/test_settings.py
```

To run a specific test:

```bash
pytest tests/test_settings.py::TestSettings::test_init_with_new_file
```

To run tests for a specific primitive operation:

```bash
# Run all featurize tests
pytest tests/test_primitives.py::TestFeaturize

# Run all common utility tests
pytest tests/test_primitives.py::TestPrimitivesCommonUtils

# Run a specific primitive test
pytest tests/test_primitives.py::TestFeaturize::test_featurize_success
```

## Test Categories

### Unit Tests

Unit tests mock external dependencies and test individual components in isolation. These tests:
- Use the `responses` library to mock HTTP requests
- Create temporary files for testing file operations
- Don't require external services

### Integration Tests

Integration tests marked with `@pytest.mark.integration` require:
- A running DeepChem server (typically at http://0.0.0.0:8000)
- Access to test data files in the assets directory
- Network connectivity

To skip integration tests during regular development:

```bash
pytest -m "not integration"
```

## Test Fixtures

The `conftest.py` file provides shared fixtures:

- `temp_settings_file` - Temporary settings file for testing
- `settings_instance` - Pre-configured Settings instance
- `data_client_instance` - Data client with test settings
- `primitives_client_instance` - Primitives client with test settings
- `mock_responses` - Mock HTTP responses
- `sample_csv_file` - Temporary CSV file for testing
- `sample_sdf_file` - Temporary SDF file for testing

## Test Assets

Test files are located in the `assets/` directory:
- `zinc10.csv` - Sample CSV file for upload testing

## Configuration

Test configuration is managed in `pytest.ini`:
- Test discovery patterns
- Coverage settings
- Markers for different test types
- Output formatting

## Best Practices

1. **Isolation**: Each test should be independent and not rely on other tests
2. **Cleanup**: Use fixtures that automatically clean up temporary files
3. **Mocking**: Mock external dependencies for unit tests
4. **Descriptive Names**: Use clear, descriptive test names
5. **Documentation**: Include docstrings for complex test scenarios
6. **Organization**: Group related tests into focused test classes (as demonstrated in `test_primitives.py`)
7. **Separation of Concerns**: Keep common utilities separate from operation-specific tests

## Troubleshooting

### Import Errors

If you encounter import errors, ensure you're running tests from the correct directory:

```bash
cd /path/to/pyds
pytest
```

### Missing Dependencies

If tests fail due to missing dependencies:

```bash
pip install -r requirements.txt
```

### Integration Test Failures

If integration tests fail:
1. Ensure the DeepChem server is running
2. Check the server URL in test configuration
3. Verify test data files exist in assets/