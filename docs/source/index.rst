DeepChem Server Documentation
==============================

**What is DeepChem Server?**

DeepChem Server is a minimal cloud infrastructure for molecular machine learning that provides a FastAPI-based backend for 
managing datasets, running featurization tasks, and building machine learning models with DeepChem.

DeepChem Server simplifies molecular machine learning workflows by providing:

* **FastAPI Backend**: Modern, fast web framework with automatic API documentation
* **DeepChem Integration**: Built-in support for molecular featurization and modeling powered by DeepChem library
* **Flexible Storage**: Disk-based datastore with support for various data formats
* **Python Client**: Easy-to-use pyds library for programmatic access
* **Containerized Deployment**: First-class Docker support for easy setup and scaling

.. image:: ./assets/img/deepchem-server.png
   :align: center
   :alt: Database Schema Diagram
   :width: 800


This documentation will guide you through setting up and using DeepChem Server for your molecular machine learning projects. 
Whether you're new to the platform or looking to integrate it into existing workflows, you'll find comprehensive guides covering 
installation, server configuration, available computation primitives, and the Python client library for programmatic access.

.. toctree::
   :maxdepth: 3
   :caption: Getting Started

   get_started/prerequisites
   get_started/installation
   get_started/quick_start
   get_started/next_steps

.. toctree::
   :maxdepth: 3
   :caption: Backend Server

   backend_server/overview
   backend_server/rest_api
   backend_server/core_modules
   backend_server/utils

.. toctree::
   :maxdepth: 3
   :caption: Core Primitives

   primitives

.. toctree::
   :maxdepth: 3
   :caption: PyDS Library

   py_ds_library/getting_started
   py_ds_library/core_modules
   py_ds_library/primitives
   py_ds_library/configuration
