#!/bin/bash

# lint.sh - Run yapf, isort, flake8, and mypy on Python files in a specified directory
# Usage: ./lint.sh [directory] [options]
#
# Options:
#   --files      Specify individual Python files to lint (space-separated)
#   --fix        Apply yapf formatting and isort import sorting automatically
#   --strict     Use strict mypy checking (overrides setup.cfg)
#   --linter     Specify which linter(s) to run (comma-separated)
#                Available: yapf,isort,flake8,mypy (default: all)
#   --help       Show this help message

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
SPECIFIC_FILES=""
LINTERS="all"

# Outputs
FAILED_LINTERS=()

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
    echo "  --files      Specify individual Python files to lint (space-separated)"
    echo "  --fix        Apply yapf formatting and isort import sorting automatically"
    echo "  --strict     Use strict mypy checking (overrides setup.cfg)"
    echo "  --linter     Specify which linter(s) to run (comma-separated)"
    echo "               Available: yapf,isort,flake8,mypy (default: all)"
    echo "  --help       Show this help message"
    echo ""
    echo "Configuration:"
    echo "  The script looks for setup.cfg in the current directory or parent directory."
    echo "  If found, it uses the [yapf], [isort], [flake8], and [mypy] sections for configuration."
    echo ""
    echo "Examples:"
    echo "  $0                                          # Run all linters on current directory"
    echo "  $0 pyds                                     # Run all linters on pyds directory"
    echo "  $0 --linter yapf,flake8                     # Run only yapf and flake8"
    echo "  $0 --files file1.py --linter mypy          # Run only mypy on specific files"
    echo "  $0 --files core/file1.py file2.py --fix    # Fix formatting on specific files"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --fix)
            FIX_FORMAT=true
            shift
            ;;
        --files)
            shift
            # Collect all files until the next option or end
            while [[ $# -gt 0 && ! $1 =~ ^-- ]]; do
                if [[ "$SPECIFIC_FILES" == "" ]]; then
                    SPECIFIC_FILES="$1"
                else
                    SPECIFIC_FILES="$SPECIFIC_FILES $1"
                fi
                shift
            done
            ;;
        --linter)
            shift
            LINTERS="$1"
            # Validate linter names
            for linter in ${LINTERS//,/ }; do
                if [[ ! "$linter" =~ ^(yapf|isort|flake8|mypy|all)$ ]]; then
                    echo -e "${RED}Error: Invalid linter '$linter'${NC}"
                    echo "Available linters: yapf, isort, flake8, mypy"
                    exit 1
                fi
            done
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

# Validate specific files if provided
if [ -n "$SPECIFIC_FILES" ]; then
    PYTHON_FILES=""
    for file in $SPECIFIC_FILES; do
        if [ ! -f "$file" ]; then
            echo -e "${RED}Error: File '$file' does not exist${NC}"
            exit 1
        fi
        if [[ ! "$file" =~ \.py$ ]]; then
            echo -e "${RED}Error: '$file' is not a Python file${NC}"
            exit 1
        fi
        if [ -z "$PYTHON_FILES" ]; then
            PYTHON_FILES="$file"
        else
            PYTHON_FILES="$PYTHON_FILES $file"
        fi
    done
else
    # Check if directory exists
    if [ ! -d "$DIRECTORY" ]; then
        echo -e "${RED}Error: Directory '$DIRECTORY' does not exist${NC}"
        exit 1
    fi

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
fi

if [ -z "$PYTHON_FILES" ]; then
    echo -e "${YELLOW}No Python files found in '$DIRECTORY'${NC}"
    exit 0
fi

echo -e "${BLUE}Found Python files in '$DIRECTORY':${NC}"
echo "$PYTHON_FILES" | sed 's/^/  /'
echo ""

# Count files
FILE_COUNT=$(echo "$PYTHON_FILES" | wc -l | xargs)
echo -e "${BLUE}Total files to process: $FILE_COUNT${NC}"
echo ""

# Function to check if required tools are installed
check_tool() {
    if ! command -v "$1" &> /dev/null; then
        echo -e "${RED}Error: $1 is not installed${NC}"
        echo "Please install it with: pip install $1"
        exit 1
    fi
}

# Function to check if a linter should run
should_run_linter() {
    local linter=$1
    if [ "$LINTERS" = "all" ]; then
        return 0
    fi
    if [[ "$LINTERS" =~ (^|,)$linter(,|$) ]]; then
        return 0
    fi
    return 1
}

# Check if required tools are installed (only for selected linters)
echo -e "${BLUE}Checking required tools...${NC}\n"
if should_run_linter "yapf"; then
    check_tool yapf
fi
if should_run_linter "isort"; then
    check_tool isort
fi
if should_run_linter "flake8"; then
    check_tool flake8
fi
if should_run_linter "mypy"; then
    check_tool mypy
fi

# Run yapf (formatter)
if should_run_linter "yapf"; then
    echo -e "${BLUE}Running yapf (code formatter)...${NC}"
    if [ "$FIX_FORMAT" = true ]; then
        echo -e "${YELLOW}Auto-fixing formatting issues...${NC}"
        echo "$PYTHON_FILES" | xargs yapf --in-place
        echo -e "${GREEN}✓ Formatting applied${NC}"
    else
        YAPF_DIFF=$(echo "$PYTHON_FILES" | xargs yapf --diff || true)
        if [ -n "$YAPF_DIFF" ]; then
            echo -e "\n$YAPF_DIFF"
            echo -e "${RED}✗ Yapf formatting issues found:${NC}"
            echo -e "Run with --fix to automatically fix formatting issues"
            YAPF_FAILED=true
            FAILED_LINTERS+=("yapf")
        else
            echo -e "${GREEN}✓ All files are properly formatted${NC}"
            YAPF_FAILED=false
        fi
    fi
    echo ""
else
    YAPF_FAILED=false
fi

# Run isort (import sorter)
if should_run_linter "isort"; then
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
            FAILED_LINTERS+=("isort")
        fi
    fi
    echo ""
else
    ISORT_FAILED=false
fi

# Run flake8 (linter)
if should_run_linter "flake8"; then
    echo -e "${BLUE}Running flake8 (linter)...${NC}"
    if echo "$PYTHON_FILES" | xargs flake8 --count; then
        echo -e "${GREEN}✓ No linting issues found${NC}"
        FLAKE8_FAILED=false
    else
        echo -e "${RED}✗ Linting issues found${NC}"
        FLAKE8_FAILED=true
        FAILED_LINTERS+=("flake8")
    fi
    echo ""
else
    FLAKE8_FAILED=false
fi

# Run mypy (type checker)
if should_run_linter "mypy"; then
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
        FAILED_LINTERS+=("mypy")
    fi
    echo ""
else
    MYPY_FAILED=false
fi

# Summary
echo -e "${BLUE}===== SUMMARY =====${NC}"
if [ -n "$SPECIFIC_FILES" ]; then
    echo -e "Files processed: $FILE_COUNT specific file(s)"
else
    echo -e "Directory: $DIRECTORY"
    echo -e "Files processed: $FILE_COUNT"
fi
echo -e "Selected linters: $LINTERS"
echo ""

if [ "$FIX_FORMAT" = true ]; then
    if should_run_linter "yapf"; then
        echo -e "yapf (formatter): ${GREEN}APPLIED${NC}"
    fi
    if should_run_linter "isort"; then
        echo -e "isort (import sorter): ${GREEN}APPLIED${NC}"
    fi
else
    if should_run_linter "yapf"; then
        if [ "$YAPF_FAILED" = true ]; then
            echo -e "yapf (formatter): ${RED}FAILED${NC}"
        else
            echo -e "yapf (formatter): ${GREEN}PASSED${NC}"
        fi
    fi

    if should_run_linter "isort"; then
        if [ "$ISORT_FAILED" = true ]; then
            echo -e "isort (import sorter): ${RED}FAILED${NC}"
        else
            echo -e "isort (import sorter): ${GREEN}PASSED${NC}"
        fi
    fi
fi

if should_run_linter "flake8"; then
    if [ "$FLAKE8_FAILED" = true ]; then
        echo -e "flake8 (linter): ${RED}FAILED${NC}"
    else
        echo -e "flake8 (linter): ${GREEN}PASSED${NC}"
    fi
fi

if should_run_linter "mypy"; then
    if [ "$MYPY_FAILED" = true ]; then
        echo -e "mypy (type checker): ${RED}FAILED${NC}"
    else
        echo -e "mypy (type checker): ${GREEN}PASSED${NC}"
    fi
fi

# Exit with error code if any tool failed (unless we're just fixing formatting)
if [ "$FIX_FORMAT" = false ]; then
    if [ "$YAPF_FAILED" = true ] || [ "$ISORT_FAILED" = true ] || [ "$FLAKE8_FAILED" = true ] || [ "$MYPY_FAILED" = true ]; then
        echo -e "\n${RED}One or more linters failed: ${FAILED_LINTERS[*]}${NC}"
        echo "Please fix the issues and try again."
        exit 1
    fi
fi

echo ""
echo -e "${GREEN}All checks passed!${NC}"
