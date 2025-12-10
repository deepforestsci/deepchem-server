# DeepChem Core Workflow

Workflow execution engine for DeepChem Server that orchestrates primitives into compute workflows.

## Features

- **ComputeWorkflow**: Execute computational primitives based on program configuration
- **Program Map**: Registry of available computational programs (featurize, train, evaluate, etc.)

## Installation

### Development (UV Workspace)

When working within the monorepo, use UV workspace:

```bash
# From repository root
uv sync
```

### Standalone Installation

```bash
pip install deepchem-core-workflow
```

This will automatically install `deepchem-core-common` as a dependency.

## Dependencies

This package depends on:
- `deepchem-core-common`: Common utilities (address, cards, datastore, config)
- `deepchem-core-primitives`: Computational primitives (should be installed via conda)

## Usage

```python
from deepchem_server.core.workflow import ComputeWorkflow, program_map

# Define a program configuration
program = {
    'program_name': 'featurize',
    'dataset_address': 'deepchem://profile/project/data.csv',
    'featurizer': 'ecfp',
    'output': 'featurized_data',
    'dataset_column': 'smiles',
    'label_column': 'y',
}

# Create and execute workflow
workflow = ComputeWorkflow(program)
result = workflow.execute()
print(f"Output: {result}")
```

## Available Programs

The workflow engine supports the following programs:

| Program Name | Description | Module |
|-------------|-------------|--------|
| `featurize` | Molecular featurization | `primitives.feat` |
| `train` | Model training | `primitives.train` |
| `evaluate` | Model evaluation | `primitives.evaluator` |
| `infer` | Model inference | `primitives.inference` |
| `train_valid_test_split` | Dataset splitting | `primitives.splitter` |
| `generate_pose` | Molecular docking | `primitives.docking` |
| `relative_binding_free_energy` | RBFE calculations | `primitives.fep.rbfe` |
| `collate_rbfe_results` | Collate FEP results | `primitives.fep.rbfe` |

## Development

```bash
# Install with dev dependencies
uv pip install -e ".[dev]"

# Run tests
pytest tests/
```

## Package Structure

```
workflow/
├── __init__.py
├── pyproject.toml
├── README.md
├── compute.py          # ComputeWorkflow and program_map
└── tests/
    ├── __init__.py
    └── test_compute.py
```

