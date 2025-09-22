#!/bin/bash

# Default values
export DATA_DIR=./local_datastore
export DOCKERFILE=Dockerfile

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -f|--dockerfile)
      export DOCKERFILE="$2"
      shift # past argument
      shift # past value
      ;;
    -d|--data-dir)
      export DATA_DIR="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      echo "Usage: $0 [OPTIONS]"
      echo "Options:"
      echo "  -f, --dockerfile DOCKERFILE    Specify dockerfile to use (default: Dockerfile)"
      echo "  -d, --data-dir DIRECTORY       Specify data directory (default: ./local_datastore)"
      echo "  -h, --help                     Show this help message"
      echo ""
      echo "Examples:"
      echo "  $0                             # Use default Dockerfile"
      echo "  $0 -f Dockerfile.gpu          # Use GPU dockerfile"
      echo "  $0 -f Dockerfile.gpu -d /path/to/data  # Use GPU dockerfile with custom data dir"
      exit 0
      ;;
    *)
      echo "Unknown option $1"
      echo "Use -h or --help for usage information"
      exit 1
      ;;
  esac
done

echo "Building DeepChem server with dockerfile: $DOCKERFILE"
echo "Data directory: $DATA_DIR"

# Build Deepchem server image
docker compose build

# Run deepchem server
docker compose up -d