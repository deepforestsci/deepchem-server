"""
DeepChem Server Client Package

A Python client for interacting with the DeepChem server API.
Provides Primitives client for computation tasks and Data client for data operations.
"""

__version__ = "0.1.0"

from .base import BaseClient
from .data import Data
from .primitives.base import Primitive
from .primitives.evaluate import Evaluate
from .primitives.featurize import Featurize
from .primitives.docking import Docking
from .primitives.infer import Infer
from .primitives.splitter import TVTSplit
from .primitives.train import Train
from .settings import Settings


__all__ = [
    "Settings",
    "Data",
    "BaseClient",
    "Primitive",
    "Featurize",
    "Docking",
    "Train",
    "Evaluate",
    "Infer",
    "TVTSplit",
]
