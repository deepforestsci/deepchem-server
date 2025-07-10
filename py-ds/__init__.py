"""
DeepChem Data Science Client Package

A Python client for interacting with the DeepChem server API.
Provides Primitives client for computation tasks and Data client for data operations.
"""

__version__ = "0.1.0"

from .settings import Settings
from .primitives import Primitives
from .data import Data
from .base import BaseClient

__all__ = ["Settings", "Primitives", "Data", "BaseClient"]
