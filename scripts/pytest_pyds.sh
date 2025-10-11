#!/bin/bash

set -e  # Exit on any error

# Project paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
PYDS_DIR="$PROJECT_ROOT/pyds"

# Check if pyds directory exists
if [[ ! -d "$PYDS_DIR" ]]; then
    echo "Error: pyds directory not found at $PYDS_DIR"
    echo "Please run this script from the project root"
    exit 1
fi

cd "$PYDS_DIR"

# Check if pytest is available
if ! command -v python -m pytest &> /dev/null; then
    echo "Error: pytest not found. Install with: pip install pytest"
    exit 1
fi

# Check if pytest-cov is available
if ! python -c "import pytest_cov" &> /dev/null; then
    echo "Error: pytest-cov not found. Install with: pip install pytest-cov"
    exit 1
fi

# Run tests with coverage
echo "Running pyds tests with coverage..."
python -m pytest --cov --cov-config=.coveragerc --cov-report=term tests/
