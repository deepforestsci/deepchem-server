# DeepChem Core Primitives

Computational primitives for DeepChem Server including:

- **Featurization** (`feat.py`): Molecular featurization using ECFP, GraphConv, Weave, etc.
- **Training** (`train.py`): Model training with various DeepChem models
- **Evaluation** (`evaluator.py`): Model evaluation and metrics
- **Inference** (`inference.py`): Running predictions on new data
- **Splitting** (`splitter.py`): Train/validation/test dataset splitting
- **Docking** (`docking.py`): Molecular docking using VINA
- **Merge** (`merge.py`): Dataset merging utilities
- **FEP** (`fep/`): Free Energy Perturbation calculations (RBFE)

## Installation

This package requires conda/mamba due to dependencies on conda-forge packages
(pdbfixer, vina, openfe, openmm) that are not available via pip.

### Using Mamba (Recommended)

```bash
# Create and activate environment
mamba env create -f environment.yml
conda activate deepchem-core-primitives

# Install package in development mode
pip install -e .
```

### Using Conda

```bash
# Create and activate environment
conda env create -f environment.yml
conda activate deepchem-core-primitives

# Install package in development mode
pip install -e .
```

## Dependencies

### Conda-only Dependencies
These packages are only available via conda-forge:
- `pdbfixer` - Protein structure preparation
- `vina` - Molecular docking
- `openfe` - Free energy calculations
- `openmm` - Molecular dynamics

### Pip Dependencies
- `deepchem` - Core ML framework for chemistry
- `rdkit` - Cheminformatics toolkit
- `torch` - PyTorch for deep learning
- `transformers` - Hugging Face transformers
- `pandas`, `numpy`, `scikit-learn` - Data science stack

## Usage

```python
from deepchem_server.core.primitives import (
    featurize,
    train,
    model_evaluator,
    infer,
    train_valid_test_split,
    generate_pose,
)

# Featurize a dataset
feat_address = featurize(
    dataset_address="deepchem://profile/project/data.csv",
    featurizer="ecfp",
    output="featurized_data",
    dataset_column="smiles",
)

# Train a model
model_address = train(
    model_type="gcn",
    dataset_address=feat_address,
    model_name="my_model",
)

# Run inference
predictions_address = infer(
    model_address=model_address,
    data_address=feat_address,
    output="predictions.csv",
)
```

## Development

```bash
# Install with dev dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/
```

## Package Structure

```
primitives/
├── __init__.py
├── environment.yml      # Conda environment definition
├── setup.py            # Package setup (for pip install -e .)
├── README.md
├── feat.py             # Featurization
├── train.py            # Model training
├── evaluator.py        # Model evaluation
├── inference.py        # Inference/prediction
├── splitter.py         # Dataset splitting
├── docking.py          # Molecular docking
├── merge.py            # Dataset merging
├── model_mappings.py   # Model type mappings
├── model_config_mapper.py  # Model configuration
├── fep/                # Free Energy Perturbation
│   ├── __init__.py
│   └── rbfe/           # Relative Binding Free Energy
│       ├── run_rbfe.py
│       ├── system_setup.py
│       ├── collate_rbfe_results.py
│       └── utils/
└── tests/              # Package tests
    ├── conftest.py
    ├── test_feat.py
    ├── test_train.py
    └── ...
```

