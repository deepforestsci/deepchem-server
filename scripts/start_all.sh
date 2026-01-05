#!/bin/bash
# Start All Services
# This script starts Redis, Datastore, Gateway, and Workers in the correct order
# Use this for local development

set -e

cd "$(dirname "$0")/.."

# Environment variables
export REDIS_URL=${REDIS_URL:-redis://localhost:6379/0}
export JOBS_DATABASE_URL=${JOBS_DATABASE_URL:-sqlite:///./data/jobs.db}
export DATASTORE_URL=${DATASTORE_URL:-http://localhost:8081}
export DATASTORE_API_KEY=${DATASTORE_API_KEY:-dev-api-key}
export DATASTORE_BASE_DIR=${DATASTORE_BASE_DIR:-./data/datastore}
export DATASTORE_PORT=${DATASTORE_PORT:-8081}
export DATADIR=${DATADIR:-./data}

# macOS fork safety
export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES

# Create directories
mkdir -p "$DATADIR"
mkdir -p "$DATASTORE_BASE_DIR"

# Number of workers (default: 2)
NUM_WORKERS=${NUM_WORKERS:-2}

echo "=========================================="
echo "Starting DeepChem Server Stack"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  Redis:     $REDIS_URL"
echo "  Datastore: $DATASTORE_URL"
echo "  Gateway:   http://localhost:8000"
echo "  Jobs DB:   $JOBS_DATABASE_URL"
echo "  Workers:   $NUM_WORKERS"
echo ""

# Function to cleanup background processes on exit
cleanup() {
    echo ""
    echo "Shutting down services..."
    kill $(jobs -p) 2>/dev/null || true
    wait
    echo "All services stopped."
}
trap cleanup EXIT INT TERM

# Check if Redis is running
echo "Checking Redis..."
if redis-cli ping > /dev/null 2>&1; then
    echo "  Redis is already running ✓"
else
    echo "  Starting Redis..."
    redis-server --daemonize yes
    sleep 1
    if redis-cli ping > /dev/null 2>&1; then
        echo "  Redis started ✓"
    else
        echo "  ERROR: Failed to start Redis. Install with: brew install redis"
        exit 1
    fi
fi

# Check if Datastore is running
echo "Checking Datastore..."
if curl -s "http://localhost:$DATASTORE_PORT/api/v1/healthcheck" > /dev/null 2>&1; then
    echo "  Datastore is already running ✓"
else
    echo "  Starting Datastore..."
    python -m deepchem_server.services.datastore &
    DATASTORE_PID=$!
    
    # Wait for datastore to be ready
    for i in {1..30}; do
        if curl -s "http://localhost:$DATASTORE_PORT/api/v1/healthcheck" > /dev/null 2>&1; then
            echo "  Datastore started ✓"
            break
        fi
        if [ $i -eq 30 ]; then
            echo "  ERROR: Datastore failed to start"
            exit 1
        fi
        sleep 0.5
    done
fi

# Start Gateway
echo "Starting Gateway..."
uvicorn deepchem_server.api.main:app --host 0.0.0.0 --port 8000 --reload &
GATEWAY_PID=$!
sleep 2
echo "  Gateway started ✓"

# Start Workers
echo "Starting $NUM_WORKERS worker(s)..."
for i in $(seq 1 $NUM_WORKERS); do
    python -m deepchem_server.services.jobs.run_worker --simple &
    echo "  Worker $i started ✓"
    sleep 1
done

echo ""
echo "=========================================="
echo "All services started!"
echo "=========================================="
echo ""
echo "API Documentation: http://localhost:8000/docs"
echo ""
echo "Press Ctrl+C to stop all services."
echo ""

# Wait for background processes
wait
