#!/bin/bash

# lint.sh - Run yapf, isort, flake8, and mypy on Python files in a specified directory
# Usage: ./lint.sh [directory] [options]
#
# Options:
#   --fix         Apply yapf formatting and isort import sorting automatically
#   --strict      Use strict mypy checking
#   --help        Show this help message

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
DIRECTORY="."
FIX_FORMAT=false
STRICT_MYPY=false

# Function to show usage
show_help() {
    echo "Usage: $0 [directory] [options]"
    echo ""
    echo "Run yapf, isort, flake8, and mypy on Python files in the specified directory."
    echo "The script automatically uses configuration from setup.cfg if available."
    echo ""
    echo "Arguments:"
    echo "  directory     Directory to lint (default: current directory)"
    echo ""
    echo "Options:"
    echo "  --fix         Apply yapf formatting and isort import sorting automatically"
    echo "  --strict      Use strict mypy checking (overrides setup.cfg)"
    echo "  --help        Show this help message"
    echo ""
    echo "Configuration:"
    echo "  The script looks for setup.cfg in the current directory or parent directory."
    echo "  If found, it uses the [yapf], [isort], [flake8], and [mypy] sections for configuration."
    echo ""
    echo "Examples:"
    echo "  $0                    # Lint current directory"
    echo "  $0 pyds              # Lint pyds directory"
    echo "  $0 pyds --fix        # Lint and auto-fix formatting in pyds"
    echo "  $0 . --strict        # Lint current directory with strict mypy"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --fix)
            FIX_FORMAT=true
            shift
            ;;
        --strict)
            STRICT_MYPY=true
            shift
            ;;
        --help|-h)
            show_help
            exit 0
            ;;
        -*)
            echo "Unknown option $1"
            show_help
            exit 1
            ;;
        *)
            DIRECTORY="$1"
            shift
            ;;
    esac
done

# Check if directory exists
if [ ! -d "$DIRECTORY" ]; then
    echo -e "${RED}Error: Directory '$DIRECTORY' does not exist${NC}"
    exit 1
fi

# Check if required tools are installed
check_tool() {
    if ! command -v "$1" &> /dev/null; then
        echo -e "${RED}Error: $1 is not installed${NC}"
        echo "Please install it with: pip install $1"
        exit 1
    fi
}

echo -e "${BLUE}Checking required tools...${NC}"
check_tool yapf
check_tool isort
check_tool flake8
check_tool mypy

# Check for setup.cfg file
SETUP_CFG_PATH=""
if [ -f "setup.cfg" ]; then
    SETUP_CFG_PATH="./setup.cfg"
    echo -e "${GREEN}✓ Found setup.cfg in current directory${NC}"
elif [ -f "../setup.cfg" ]; then
    SETUP_CFG_PATH="../setup.cfg"
    echo -e "${YELLOW}⚠ Found setup.cfg in parent directory${NC}"
else
    echo -e "${YELLOW}⚠ No setup.cfg found - using default configurations${NC}"
fi

if [ -n "$SETUP_CFG_PATH" ]; then
    echo -e "${BLUE}Using configuration from: $SETUP_CFG_PATH${NC}"
fi
echo ""

# Find Python files (exclude common build/cache directories)
PYTHON_FILES=$(find "$DIRECTORY" -name "*.py" -type f \
    | grep -v __pycache__ \
    | grep -v "/build/" \
    | grep -v "/dist/" \
    | grep -v "/.pytest_cache/" \
    | grep -v "/htmlcov/" \
    | grep -v "/.venv/" \
    | grep -v "/venv/" \
    | grep -v "/env/" \
    | grep -v "/.git/" \
    | sort)

if [ -z "$PYTHON_FILES" ]; then
    echo -e "${YELLOW}No Python files found in '$DIRECTORY'${NC}"
    exit 0
fi

echo -e "${BLUE}Found Python files in '$DIRECTORY':${NC}"
echo "$PYTHON_FILES" | sed 's/^/  /'
echo ""

# Count files
FILE_COUNT=$(echo "$PYTHON_FILES" | wc -l)
echo -e "${BLUE}Total files to process: $FILE_COUNT${NC}"
echo ""

# Run yapf (formatter)
echo -e "${BLUE}Running yapf (code formatter)...${NC}"
if [ "$FIX_FORMAT" = true ]; then
    echo -e "${YELLOW}Auto-fixing formatting issues...${NC}"
    echo "$PYTHON_FILES" | xargs yapf --in-place
    echo -e "${GREEN}✓ Formatting applied${NC}"
else
    if echo "$PYTHON_FILES" | xargs yapf --diff | grep -q "^---"; then
        echo -e "${RED}✗ Formatting issues found${NC}"
        echo "Run with --fix to automatically fix formatting issues"
        YAPF_FAILED=true
    else
        echo -e "${GREEN}✓ All files are properly formatted${NC}"
        YAPF_FAILED=false
    fi
fi
echo ""

# Run isort (import sorter)
echo -e "${BLUE}Running isort (import sorter)...${NC}"
if [ "$FIX_FORMAT" = true ]; then
    echo -e "${YELLOW}Auto-fixing import order...${NC}"
    echo "$PYTHON_FILES" | xargs isort
    echo -e "${GREEN}✓ Import sorting applied${NC}"
    ISORT_FAILED=false
else
    if echo "$PYTHON_FILES" | xargs isort --check-only --diff; then
        echo -e "${GREEN}✓ All imports are properly sorted${NC}"
        ISORT_FAILED=false
    else
        echo -e "${RED}✗ Import sorting issues found${NC}"
        echo "Run with --fix to automatically fix import order"
        ISORT_FAILED=true
    fi
fi
echo ""

# Run flake8 (linter)
echo -e "${BLUE}Running flake8 (linter)...${NC}"
if echo "$PYTHON_FILES" | xargs flake8; then
    echo -e "${GREEN}✓ No linting issues found${NC}"
    FLAKE8_FAILED=false
else
    echo -e "${RED}✗ Linting issues found${NC}"
    FLAKE8_FAILED=true
fi
echo ""

# Run mypy (type checker)
echo -e "${BLUE}Running mypy (type checker)...${NC}"
MYPY_ARGS=""
if [ "$STRICT_MYPY" = true ]; then
    echo -e "${YELLOW}Using strict type checking...${NC}"
    MYPY_ARGS="--strict"
fi

if echo "$PYTHON_FILES" | xargs mypy $MYPY_ARGS; then
    echo -e "${GREEN}✓ No type checking issues found${NC}"
    MYPY_FAILED=false
else
    echo -e "${RED}✗ Type checking issues found${NC}"
    MYPY_FAILED=true
fi
echo ""

# Summary
echo -e "${BLUE}=== SUMMARY ===${NC}"
echo -e "Directory: $DIRECTORY"
echo -e "Files processed: $FILE_COUNT"
echo ""

if [ "$FIX_FORMAT" = true ]; then
    echo -e "yapf (formatter): ${GREEN}Applied${NC}"
    echo -e "isort (import sorter): ${GREEN}Applied${NC}"
else
    if [ "$YAPF_FAILED" = true ]; then
        echo -e "yapf (formatter): ${RED}FAILED${NC}"
    else
        echo -e "yapf (formatter): ${GREEN}PASSED${NC}"
    fi

    if [ "$ISORT_FAILED" = true ]; then
        echo -e "isort (import sorter): ${RED}FAILED${NC}"
    else
        echo -e "isort (import sorter): ${GREEN}PASSED${NC}"
    fi
fi

if [ "$FLAKE8_FAILED" = true ]; then
    echo -e "flake8 (linter): ${RED}FAILED${NC}"
else
    echo -e "flake8 (linter): ${GREEN}PASSED${NC}"
fi

if [ "$MYPY_FAILED" = true ]; then
    echo -e "mypy (type checker): ${RED}FAILED${NC}"
else
    echo -e "mypy (type checker): ${GREEN}PASSED${NC}"
fi

# Exit with error code if any tool failed (unless we're just fixing formatting)
if [ "$FIX_FORMAT" = false ]; then
    if [ "$YAPF_FAILED" = true ] || [ "$ISORT_FAILED" = true ] || [ "$FLAKE8_FAILED" = true ] || [ "$MYPY_FAILED" = true ]; then
        echo ""
        echo -e "${RED}Some checks failed. Please fix the issues above.${NC}"
        exit 1
    fi
fi

echo ""
echo -e "${GREEN}All checks passed!${NC}"
