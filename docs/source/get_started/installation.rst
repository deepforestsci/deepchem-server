Installation
============

DeepChem Server can be installed and run using Docker (recommended) or manually from source. Docker provides the simplest setup with all dependencies pre-configured.

Docker Installation (Recommended)
---------------------------------

Docker provides containerized deployment with all dependencies managed automatically.

Prerequisites
~~~~~~~~~~~~~

* Docker installed on your system
* Docker Compose (included with Docker Desktop)

Steps
~~~~~

1. **Clone the Repository**

   .. code-block:: bash

      git clone <repository-url>
      cd deepchem-server

2. **Start the Server**

   Use the provided helper script for easy setup:

   .. code-block:: bash

      ./docker.sh

   **Available Options:**

   .. code-block:: bash

      ./docker.sh

   **Examples:**

   .. code-block:: bash

      # Use default CPU setup
      ./docker.sh

      # Use GPU setup
      ./docker.sh -f Dockerfile.gpu

3. **Verify Installation**

   Check that the server is running:

   .. code-block:: bash

      curl http://localhost:8000/healthcheck

   Expected response: ``{"status": "ok"}``

   Access the interactiveAPI documentation at: http://localhost:8000/docs

Available Dockerfiles
~~~~~~~~~~~~~~~~~~~~~~

* **Dockerfile**: CPU-based setup using micromamba base image
* **Dockerfile.gpu**: GPU-accelerated setup using NVIDIA PyTorch base image

**Dockerfile Features:**

* Uses micromamba for fast package management
* Includes all DeepChem dependencies
* Optimized for CPU-based molecular machine learning
* Health check endpoint configured

**Dockerfile.gpu Features:**

* NVIDIA PyTorch base image with CUDA support
* GPU-accelerated DeepChem operations
* Compatible with NVIDIA Docker runtime
* Optimized for large-scale molecular modeling

Manual Installation
-------------------

For development or custom deployments, install DeepChem Server manually using micromamba.

Prerequisites
~~~~~~~~~~~~~

* micromamba package manager (lightweight conda replacement)
* Python 3.11 (specified in environment file)
* Git for cloning the repository

Steps
~~~~~

1. **Clone the Repository**

   .. code-block:: bash

      git clone <repository-url>
      cd deepchem-server

2. **Create Environment from YAML**

   Use the provided environment file to create a consistent environment:

   .. code-block:: bash

      micromamba env create -f deepchem_server/environments/core_environment.yml

3. **Activate Environment**

   .. code-block:: bash

      micromamba activate deepchem-server-env

4. **Start the Server**

   For development with auto-reload:

   .. code-block:: bash

      cd deepchem_server
      uvicorn main:app --host 0.0.0.0 --port 8000 --reload

   Or use the development script:

   .. code-block:: bash

      ./start-dev-server.sh

Configuration
-------------

**Environment Management**

* **Environment File**: ``deepchem_server/environments/core_environment.yml``
* **Python Version**: 3.11
* **Package Manager**: micromamba for dependency management

**Environment Variables**

* ``DATADIR``: Data directory path (default: ``/opt/deepchem_server_app/data`` in Docker)
* ``PYTHONPATH``: Python path for module imports

**Port Configuration**

Default port is 8000. To change:

.. code-block:: bash

   # Docker
   docker run -p 8080:8000 deepchem-server

   # Manual
   uvicorn main:app --port 8080

Verification
------------

After installation, verify the setup:

1. **Health Check**

   .. code-block:: bash

      curl http://localhost:8000/healthcheck

2. **API Documentation**

   Visit http://localhost:8000/docs for interactive API documentation

3. **Test Upload**

   .. code-block:: bash

      cd pyds/tests
      python test_upload_featurize.py

Getting Help
~~~~~~~~~~~~

If you encounter issues:

1. Check server logs for error messages
2. Verify micromamba environment is correctly created and activated
3. Ensure system meets minimum requirements
4. Verify Python version is 3.11 as specified in environment file
5. Consult the FastAPI documentation for deployment questions
6. Check the DeepChem documentation for primitive operations