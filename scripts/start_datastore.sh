#!/bin/bash
# Start Datastore Service
# Provides file storage API for data and models

set -e

cd "$(dirname "$0")/.."

# Environment variables
export DATASTORE_BASE_DIR=./data/datastore
export DATASTORE_API_KEY=${DATASTORE_API_KEY:-dev-api-key}
export DATASTORE_PORT=${DATASTORE_PORT:-8081}

# Create data directory
mkdir -p "$DATASTORE_BASE_DIR"

echo "Starting Datastore Service..."
echo "  Base Dir: $DATASTORE_BASE_DIR"
echo "  Port: $DATASTORE_PORT"
echo "  URL: http://localhost:$DATASTORE_PORT"

python -m deepchem_server.services.datastore
