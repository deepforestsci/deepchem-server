# flake8: noqa
"""
DeepChem Server Core Package

This package contains three sub-packages:
- common: Common utilities (address, cards, config, datastore)
- primitives: Computational primitives (feat, train, inference, docking, fep)
- workflow: Workflow execution engine

Each sub-package can be installed independently:
- deepchem-core-common: Uses UV (pyproject.toml)
- deepchem-core-primitives: Uses conda/mamba (environment.yml + setup.py)
- deepchem-core-workflow: Uses UV (pyproject.toml)

Note: This __init__.py uses lazy imports to avoid loading heavy dependencies
(like torch, deepchem models) until they are actually needed.
"""

from typing import TYPE_CHECKING

# Always available - common sub-package (lightweight)
from deepchem_server.core.common import (
    Card,
    config,
    DataCard,
    DataStore,
    DeepchemAddress,
    DiskDataStore,
    get_datastore,
    log_progress,
    ModelCard,
    parse_boolean_none_values_from_kwargs,
    parse_dict_with_datatypes,
    run_job,
    set_datastore,
)


# Lazy imports for primitives and workflow (heavy dependencies)
# These are only loaded when explicitly accessed
def __getattr__(name: str):
    """Lazy import for primitives and workflow modules."""
    # Primitives
    if name == "generate_pose":
        from deepchem_server.core.primitives import generate_pose

        return generate_pose
    elif name == "model_evaluator":
        from deepchem_server.core.primitives import model_evaluator

        return model_evaluator
    elif name == "featurize":
        from deepchem_server.core.primitives import featurize

        return featurize
    elif name == "featurizer_map":
        from deepchem_server.core.primitives import featurizer_map

        return featurizer_map
    elif name == "infer":
        from deepchem_server.core.primitives import infer

        return infer
    elif name == "merge":
        from deepchem_server.core.primitives import merge

        return merge
    elif name == "train_valid_test_split":
        from deepchem_server.core.primitives import train_valid_test_split

        return train_valid_test_split
    elif name == "train":
        from deepchem_server.core.primitives import train

        return train
    elif name == "model_mappings":
        from deepchem_server.core.common import model_mappings

        return model_mappings
    # Workflow
    elif name == "ComputeWorkflow":
        from deepchem_server.core.workflow import ComputeWorkflow

        return ComputeWorkflow
    elif name == "program_map":
        from deepchem_server.core.workflow import program_map

        return program_map

    raise AttributeError(f"module 'deepchem_server.core' has no attribute '{name}'")
