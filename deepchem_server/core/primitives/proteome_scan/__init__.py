"""
ProteomeScan primitives for deepchem-server.

Primitives used to run the proteome scan workflow in deepchem-server.
"""

from deepchem_server.core.primitives.proteome_scan.filter_promiscuous_targets import (
    filter_promiscuous_targets,
)


__all__ = [
    "filter_promiscuous_targets",
]
