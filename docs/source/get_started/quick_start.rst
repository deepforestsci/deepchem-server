Quick Start
===========

This guide will help you run your first molecular machine learning workflow with DeepChem Server.

Before You Begin
----------------

Ensure you have:

1. DeepChem Server running (see :doc:`installation`)
2. A sample dataset (CSV file with SMILES strings or molecular data)
3. Web browser for accessing the interactive API documentation

Interactive API Documentation
-----------------------------

The best way to explore and test the API is through the interactive documentation:

**Swagger UI**: http://localhost:8000/docs

This provides:

* Complete endpoint documentation with request/response schemas
* Interactive request testing with real-time responses
* Parameter descriptions and validation rules
* Example requests and responses
* Schema definitions for all data models

Basic Workflow
--------------

The typical molecular machine learning workflow involves:

1. **Upload Data**: Submit your molecular dataset to the server
2. **Featurize**: Transform molecules into machine learning features
3. **Train**: Build machine learning models on featurized data
4. **Evaluate**: Assess model performance
5. **Infer**: Make predictions on new data

Available Endpoints
-------------------

**Data Management**
* ``POST /data/uploaddata``: Upload datasets to the datastore
* ``GET /data/{dataset_id}/download``: Download processed datasets

**Primitive Operations**
* ``POST /primitive/featurize``: Apply molecular featurization
* ``POST /primitive/train``: Train machine learning models
* ``POST /primitive/evaluate``: Evaluate model performance
* ``POST /primitive/infer``: Run inference on new data
* ``POST /primitive/train-valid-test-split``: Split datasets for training

**System**
* ``GET /healthcheck``: Check server health status


Python Client Library
---------------------

For programmatic access, use the pyds Python client library:

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

   # Upload and process data
   response = data_client.upload_data("dataset.csv")
   dataset_address = response['dataset_address']

   # Featurize
   response = featurize_client.run(
       dataset_address=dataset_address,
       featurizer="ECFP",
       output="featurized_data",
       dataset_column="smiles"
   )

   # Train model
   response = train_client.run(
       dataset_address=response['featurized_file_address'],
       model_type="random_forest_classifier",
       model_name="my_model"
   )

For detailed Python client documentation, see :doc:`../py_ds_library/index`.

Next Steps
----------

1. **Explore the API**: Visit http://localhost:8000/docs to test endpoints interactively
2. **Try Different Featurizers**: Experiment with various molecular representations
3. **Use the Python Client**: Explore the pyds library for programmatic workflows
4. **Build ML Pipelines**: Chain operations for complete machine learning workflows
5. **Consult Documentation**: Review detailed API reference and primitives documentation

Troubleshooting
---------------

**Server Not Responding**
   Check if the server is running:

   .. code-block:: bash

      curl http://localhost:8000/healthcheck

**Memory Issues**
   Large datasets may require additional memory. Consider using smaller batches or upgrading system resources

**Need More Information**
   Visit http://localhost:8000/docs for comprehensive API documentation and interactive testing