"""
DeepChem Server Client Package

A Python client for interacting with the DeepChem server API.
Provides Primitives client for computation tasks and Data client for data operations.
"""

__version__ = "0.1.0"

from .base import BaseClient
from .data import Data
from .primitives import Primitives
from .settings import Settings

__all__ = ["Settings", "Primitives", "Data", "BaseClient"]
