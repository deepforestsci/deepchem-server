# DeepChem Server Core

The core package contains three sub-packages that can be developed together (UV workspace) and published independently.

## Package Structure

```
core/
├── common/                 # Common utilities (UV workspace member)
│   ├── pyproject.toml     # Package definition
│   ├── address.py         # DeepChem address handling
│   ├── cards.py           # DataCard and ModelCard
│   ├── config.py          # Global configuration
│   ├── datastore.py       # DataStore implementations
│   ├── progress_logger.py # Progress logging utilities
│   ├── utils.py           # Helper functions
│   └── tests/             # Package-specific tests
│
├── primitives/            # Computational primitives (conda/mamba-managed)
│   ├── environment.yml    # Conda environment definition
│   ├── setup.py          # Package setup for pip install
│   ├── feat.py           # Featurization
│   ├── train.py          # Model training
│   ├── evaluator.py      # Model evaluation
│   ├── inference.py      # Inference/prediction
│   ├── splitter.py       # Dataset splitting
│   ├── docking.py        # Molecular docking (VINA)
│   ├── merge.py          # Dataset merging
│   ├── fep/              # Free Energy Perturbation
│   │   └── rbfe/         # Relative Binding Free Energy
│   └── tests/            # Package-specific tests
│
├── workflow/             # Workflow execution (UV workspace member)
│   ├── pyproject.toml    # Package definition
│   ├── compute.py        # ComputeWorkflow and program_map
│   └── tests/            # Package-specific tests
│
├── environments/         # Environment configurations
│   ├── core_environment.yml  # Core conda environment
│   └── gpu_environment.yml   # GPU-enabled environment
│
└── tests/               # Integration tests (shared assets)
    └── assets/          # Test data files
```

## Development Setup (UV Workspace)

The UV workspace allows developing all packages together while maintaining independent publishability.

```bash
# From repository root
uv sync
```

This installs all workspace members in editable mode with their dependencies linked.

## Individual Package Installation

### deepchem-core-common

```bash
pip install deepchem-core-common
```

### deepchem-core-primitives (Conda/Mamba)

Requires conda-forge packages not available via pip:

```bash
cd deepchem_server/core/primitives
mamba env create -f environment.yml
conda activate deepchem-core-primitives
pip install -e .
```

### deepchem-core-workflow

```bash
pip install deepchem-core-workflow
```

## Dependencies

| Package | Manager | Key Dependencies |
|---------|---------|------------------|
| `common` | UV | pandas, pillow, deepchem, numpy |
| `primitives` | Conda | rdkit, vina, openfe, openmm, deepchem, torch |
| `workflow` | UV | deepchem-core-common |

## Running Tests

```bash
# From repository root (runs all tests)
uv run pytest

# Individual package tests
cd deepchem_server/core/common && pytest tests/
cd deepchem_server/core/primitives && pytest tests/
cd deepchem_server/core/workflow && pytest tests/
```

## Publishing

Each package can be published independently to PyPI:

```bash
# Build and publish common
cd deepchem_server/core/common
uv build
uv publish

# Build and publish workflow
cd deepchem_server/core/workflow
uv build
uv publish
```

Note: `primitives` is published via `setup.py` since it uses conda:

```bash
cd deepchem_server/core/primitives
python -m build
twine upload dist/*
```

## Dependency Hierarchy

```
common (no internal dependencies)
    ↓
primitives (depends on common)
    ↓
workflow (depends on common)
```

## Test Assets

Test assets are shared and located in `core/tests/assets/`.

