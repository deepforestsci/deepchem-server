#!/bin/bash

# pytest_pyds.sh - Run pytest on the pyds package
# Usage: ./pytest_pyds.sh

cd pyds
python -m pytest tests/ -v --cov --cov-report=term-missing