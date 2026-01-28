#!/bin/bash
# Start Worker Service
# Executes DeepChem jobs from the queue
# Requires conda environment with DeepChem installed

set -e

cd "$(dirname "$0")/.."

# Environment variables
export REDIS_URL=${REDIS_URL:-redis://localhost:6379/0}
export JOBS_DATABASE_URL=${JOBS_DATABASE_URL:-sqlite:///./data/jobs.db}
export DATASTORE_URL=${DATASTORE_URL:-http://localhost:8081}
export DATASTORE_API_KEY=${DATASTORE_API_KEY:-dev-api-key}
export DATADIR=${DATADIR:-./data}

# macOS fork safety (required for DeepChem/RDKit)
export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES

# Create data directory
mkdir -p "$DATADIR"

echo "Starting Worker Service..."
echo "  Redis: $REDIS_URL"
echo "  Datastore: $DATASTORE_URL"
echo "  Jobs DB: $JOBS_DATABASE_URL"

# Use --simple on macOS to avoid fork issues with Objective-C
python -m deepchem_server.services.jobs.run_worker --simple
