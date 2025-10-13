Getting Started
===============

The pyds library is a Python client package for interacting with the DeepChem Server API. 
It provides a clean, object-oriented interface for managing settings, uploading data, and 
submitting primitive jobs for molecular machine learning workflows.

**What is pyds?**

pyds simplifies the process of working with molecular data by providing:

* **Unified API**: A consistent interface for all DeepChem Server operations
* **Settings Management**: Centralized configuration for profiles, projects, and server connections
* **Data Operations**: Easy upload and management of molecular datasets
* **ML Primitives**: Ready-to-use components for featurization, training, evaluation, and inference
* **Workflow Integration**: Seamless chaining of operations for complete ML pipelines


Installation
------------

Install the pyds library from source:

.. code-block:: bash

   cd pyds
   pip install -e .

For development with testing dependencies:

.. code-block:: bash

   pip install -e ".[dev]"

Architecture
------------

The pyds library follows a clean inheritance structure designed for modularity and code reuse:

.. code-block:: text

   BaseClient (base functionality)
   ├── Data (data operations)
   └── Primitive (abstract base for computation tasks)
       ├── Featurize (molecular featurization)
       ├── Train (model training)
       ├── Evaluate (model evaluation)
       ├── Infer (inference/predictions)
       └── TVTSplit (train-valid-test splitting)

**Key Design Principles:**

* **BaseClient**: Contains all common functionality like HTTP requests, configuration validation, and shared utilities
* **Inheritance-based**: Specific clients inherit from BaseClient, eliminating code duplication
* **Consistent Interface**: All clients provide the same base methods and configuration handling
* **Settings Management**: Centralized configuration through a Settings class with persistent storage in a JSON file

Quick Start
-----------

Basic workflow for using the pyds library:

1. **Configure Settings**: Set up profile, project, and server URL
2. **Initialize Clients**: Create Data and Primitive client instances
3. **Upload Data**: Use Data client to upload datasets
4. **Run Primitives**: Use primitive classes for computation tasks

.. code-block:: python

   from pyds import Settings, Data, Featurize, Train

   # Configure settings
   settings = Settings()
   settings.set_profile("my_profile")
   settings.set_project("my_project")

   # Initialize clients
   data_client = Data(settings)
   featurize_client = Featurize(settings)
   train_client = Train(settings)

   # Upload data
   response = data_client.upload_data("data.csv", description="My dataset")
   dataset_address = response['dataset_address']

   # Featurize data
   response = featurize_client.run(
       dataset_address=dataset_address,
       featurizer="ECFP",
       output="featurized_data",
       dataset_column="smiles"
   )
   featurized_address = response['featurized_file_address']

   # Train model
   response = train_client.run(
       dataset_address=featurized_address,
       model_type="random_forest_classifier",
       model_name="my_model"
   )
