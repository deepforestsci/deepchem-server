# flake8: noqa
"""
DeepChem Server Core Package.

This package contains:
- common: Shared utilities (cards, address, config) - NO ML DEPENDENCIES
- primitives: ML operations (train, evaluate, featurize) - REQUIRES DEEPCHEM

Note: primitives are NOT imported here to avoid loading deepchem at import time.
Import them explicitly where needed:
    from deepchem_server.core.primitives import train, evaluate, etc.
"""

from deepchem_server.core.common import config
from deepchem_server.core.common import cards
from deepchem_server.core.common.cards import Card, DataCard, ModelCard
from deepchem_server.core.common.address import DeepchemAddress
from deepchem_server.core.common.config import set_datastore, get_datastore
from deepchem_server.core.common.progress_logger import log_progress
from deepchem_server.core.common.model_mappings import model_address_map

# NOTE: primitives are NOT imported here to keep gateway lightweight
# Workers should import primitives explicitly:
#   from deepchem_server.core.primitives import evaluator, splitter, train, etc.
