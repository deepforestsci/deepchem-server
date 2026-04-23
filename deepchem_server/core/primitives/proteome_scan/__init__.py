"""
ProteomeScan primitives for deepchem-server.

Primitives used to run the proteome scan workflow in deepchem-server.
"""

from deepchem_server.core.primitives.proteome_scan.ligand_prep import ligand_prep


__all__ = [
    "ligand_prep",
]
