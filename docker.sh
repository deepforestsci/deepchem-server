#!/bin/bash

# Default values
export DOCKERFILE=Dockerfile
SHOW_LOGS=false
FOLLOW_LOGS=false
TAIL_LOGS=100
TAIL_EXPLICIT=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'


log_info() {
    echo -e "${BLUE}$1${NC}"
}
log_error() {
    echo -e "${RED}$1${NC}"
}
log_success() {
    echo -e "${GREEN}$1${NC}"
}
log_warning() {
    echo -e "${YELLOW}$1${NC}"
}


docker_compose_cmd() {
    if command -v docker-compose &> /dev/null; then
        echo "docker-compose"
    elif docker compose version &> /dev/null; then
        echo "docker compose"
    else
        return 1
    fi
}

DOCKER_COMPOSE_CMD=$(docker_compose_cmd)
if [ $? -ne 0 ]; then
    log_error "Error: docker compose not found"
    exit 1
fi

# Function to display docker compose logs
show_docker_logs() {
    local tail_count=${1:-$TAIL_LOGS}
    local follow=${2:-false}
    local show_all=${3:-false}
    
    log_info "=== DOCKER COMPOSE LOGS ==="
    echo "Showing logs from deepchem-server service..."
    
    # Build the command
    local cmd="$DOCKER_COMPOSE_CMD logs --timestamps"
    
    if [ "$follow" = true ]; then
        cmd="$cmd --follow"
        log_info "Following logs (press Ctrl+C to stop)..."
    elif [ "$show_all" = true ]; then
        log_info "Showing all available logs..."
    else
        cmd="$cmd --tail=$tail_count"
        log_info "Showing last $tail_count lines..."
    fi
    
    cmd="$cmd deepchem-server"
    
    log_info "Running: $cmd"
    $cmd || log_error "Failed to retrieve logs"
    
    log_info "=== END DOCKER LOGS ==="
}

# Function to check service status
check_service_status() {
    if $DOCKER_COMPOSE_CMD ps deepchem-server | grep -q "Up"; then
        log_success "✓ deepchem-server is running"
        return 0
    else
        log_warning "⚠ deepchem-server is not running"
        return 1
    fi
}

healthcheck_status() {

    max_retries=3
    retry_delay=5
    retries=1

    while [ $retries -le $max_retries ]; do
        if curl -f http://localhost:8000/healthcheck &> /dev/null; then
            log_success "✓ deepchem-server is healthy"
            return 0
        else
            log_warning "Healthcheck failed. Retrying... [$retries/$max_retries]"
            retries=$((retries + 1))
            sleep $retry_delay
        fi
    done
    return 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -f|--dockerfile)
      export DOCKERFILE="$2"
      shift
      shift
      ;;
    -l|--logs)
      SHOW_LOGS=true
      shift
      ;;
    -lf|--logs-follow)
      FOLLOW_LOGS=true
      shift
      ;;
    -lt|--logs-tail)
      TAIL_LOGS="$2"
      TAIL_EXPLICIT=true
      shift
      shift
      ;;
    -st|--stop)
      $DOCKER_COMPOSE_CMD down
      exit $?
      ;;
    -s|--status)
      check_service_status
      healthcheck_status
      exit $?
      ;;
    -h|--help)
      echo "Usage: $0 [OPTIONS]"
      echo "Options:"
      echo "  -f, --dockerfile DOCKERFILE    Specify dockerfile to use (default: Dockerfile)"
      echo "  -l, --logs                     Show logs (all available logs by default)"
      echo "  -lf, --logs-follow             Follow logs in real-time after starting"
      echo "  -lt, --logs-tail NUMBER        Number of log lines to show (default: $TAIL_LOGS)"
      echo "  -st, --stop                     Stop the deepchem server"
      echo "  -s, --status                   Check service status and exit"
      echo "  -h, --help                     Show this help message"
      echo ""
      echo "Examples:"
      echo "  $0                             # Use default Dockerfile"
      echo "  $0 -f Dockerfile.gpu           # Use GPU dockerfile"
      echo "  $0 -l                          # Show all logs"
      echo "  $0 -lf                         # Start and follow logs"
      echo "  $0 -lt 50                      # Show last 50 log lines"
      echo "  $0 -st                         # Stop the deepchem server"
      echo "  $0 -s                          # Check service status"
      echo "  $0 -f Dockerfile.gpu           # Use GPU dockerfile"
      exit 0
      ;;
    *)
      echo "Unknown option $1"
      echo "Use -h or --help for usage information"
      exit 1
      ;;
  esac
done

# Determine if we need to build and start the service
SKIP_BUILD=false

if [ "$SHOW_LOGS" = true ] && [ "$FOLLOW_LOGS" = false ]; then
    SKIP_BUILD=true
    log_info "Logs-only mode: skipping build and start"
fi

if [ "$FOLLOW_LOGS" = true ]; then
    SKIP_BUILD=true
    log_info "Follow logs mode: skipping build and start"
fi

if [ "$SKIP_BUILD" = false ]; then
    log_info "Building DeepChem server with dockerfile: $DOCKERFILE"

    log_info "Building docker image..."
    $DOCKER_COMPOSE_CMD build

    if [ $? -ne 0 ]; then
        log_error "✗ Build failed"
        exit 1
    fi

    log_success "✓ Build completed"

    log_info "Starting deepchem server..."
    $DOCKER_COMPOSE_CMD up -d

    if [ $? -ne 0 ]; then
        log_error "✗ Failed to start deepchem server"
        exit 1
    fi

    log_success "✓ Deepchem server started"

    sleep 2
    check_service_status
fi

if [ "$FOLLOW_LOGS" = true ]; then
    show_docker_logs $TAIL_LOGS true false
elif [ "$SHOW_LOGS" = true ]; then
    if [ "$TAIL_EXPLICIT" = false ]; then
        show_docker_logs $TAIL_LOGS false true
    else
        show_docker_logs $TAIL_LOGS false false
    fi
fi