#!/bin/bash
# Start Redis server for job queue
# Run this first before other services

set -e

echo "Starting Redis server..."
redis-server

# If redis-server is not installed:
# brew install redis
