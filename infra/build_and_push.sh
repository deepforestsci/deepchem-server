#!/usr/bin/env bash
# Build and push CPU and GPU compute images to ECR.
# Run from the repo root: bash infra/build_and_push.sh
set -euo pipefail

REGION=${AWS_DEFAULT_REGION:-us-east-1}
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
ECR_BASE="${ACCOUNT}.dkr.ecr.${REGION}.amazonaws.com"
REPO="${ECR_BASE}/${PROJECT_NAME:-deepchem-server}-compute"

echo "Authenticating with ECR..."
aws ecr get-login-password --region "$REGION" \
  | docker login --username AWS --password-stdin "$ECR_BASE"

echo "Building CPU image..."
docker buildx build \
  --platform linux/amd64 \
  -t "${REPO}:cpu" \
  -f infra/Dockerfile.cpu \
  . --push

echo "Building GPU image..."
docker buildx build \
  --platform linux/amd64 \
  -t "${REPO}:gpu" \
  -f infra/Dockerfile.gpu \
  . --push

echo ""
echo "Done. Images pushed:"
echo "  ${REPO}:cpu"
echo "  ${REPO}:gpu"
