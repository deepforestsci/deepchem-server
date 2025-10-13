Core Modules
============

This section provides detailed API documentation for the core modules of DeepChem Server, automatically generated from the source code docstrings.

Datastore Operations
--------------------

The datastore module handles persistent storage of datasets and models.

.. automodule:: deepchem_server.core.datastore
   :members:
   :undoc-members:
   :show-inheritance:

Address Management
------------------

The address module provides URI-like addressing for resources.

.. automodule:: deepchem_server.core.address
   :members:
   :undoc-members:
   :show-inheritance:

Cards and Metadata
------------------

The cards module defines metadata structures for datasets and models.

.. automodule:: deepchem_server.core.cards
   :members:
   :undoc-members:
   :show-inheritance:

Configuration
-------------

The config module manages server configuration settings.

.. automodule:: deepchem_server.core.config
   :members:
   :undoc-members:
   :show-inheritance:

Featurization
-------------

The feat module provides molecular featurization capabilities.

.. automodule:: deepchem_server.core.feat
   :members:
   :undoc-members:
   :show-inheritance:

Compute Operations
------------------

The compute module handles computational tasks and job execution.

.. automodule:: deepchem_server.core.compute
   :members:
   :undoc-members:
   :show-inheritance:

Molecular Docking
------------------

The docking module provides molecular docking capabilities using AutoDock VINA.

**Key Features:**
- Generates protein-ligand binding poses using AutoDock VINA
- Supports both PDB and PDBQT output formats
- Automatically splits PDBQT files for multiple binding modes
- Returns DeepChem addresses to all generated files

.. automodule:: deepchem_server.core.docking
   :members:
   :undoc-members:
   :show-inheritance: 