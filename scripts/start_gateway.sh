#!/bin/bash
# Start Gateway Service (API Server)
# Provides REST API for DeepChem operations

set -e

cd "$(dirname "$0")/.."

# Environment variables
export REDIS_URL=${REDIS_URL:-redis://localhost:6379/0}
export JOBS_DATABASE_URL=${JOBS_DATABASE_URL:-sqlite:///./data/jobs.db}
export DATASTORE_URL=${DATASTORE_URL:-http://localhost:8081}
export DATASTORE_API_KEY=${DATASTORE_API_KEY:-dev-api-key}
export DATADIR=${DATADIR:-./data}

# Create data directory
mkdir -p "$DATADIR"

echo "Starting Gateway Service..."
echo "  Redis: $REDIS_URL"
echo "  Datastore: $DATASTORE_URL"
echo "  Jobs DB: $JOBS_DATABASE_URL"
echo "  URL: http://localhost:8000"

uvicorn deepchem_server.api.main:app --host 0.0.0.0 --port 8000 --reload
