Core Primitives
===============

This section documents the core primitives that form the backbone of DeepChem Server's machine learning workflow. These primitives provide the essential functionality for molecular machine learning pipelines: data featurization, model training, inference, and evaluation.

Overview
--------

DeepChem Server provides five main primitives that work together to create end-to-end machine learning workflows:

* **Featurize**: Transform raw molecular data into machine learning features
* **Train**: Build and train machine learning models on featurized datasets
* **Inference**: Run predictions on new data using trained models
* **Evaluation**: Assess model performance using various metrics
* **Docking**: Perform molecular docking to predict protein-ligand binding poses

These primitives are designed to work seamlessly together while also being usable independently for specific tasks.

Featurization
-------------

The featurization primitive transforms raw molecular data (like SMILES strings or SDF files) into numerical features that can be used for machine learning.

.. autofunction:: deepchem_server.core.feat.featurize
   :no-index:

Supporting Functions
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: deepchem_server.core.feat.split_dataset
   :no-index:

.. autofunction:: deepchem_server.core.feat.featurize_part
   :no-index:

.. autofunction:: deepchem_server.core.feat.featurize_multi_core
   :no-index:

Available Featurizers
~~~~~~~~~~~~~~~~~~~~~

The featurization primitive supports the following DeepChem featurizers:

* **ecfp**: Extended Connectivity Fingerprints (Circular Fingerprints) - Compatible with scikit-learn models
* **graphconv**: Graph Convolution Featurizer - Compatible with scikit-learn models
* **weave**: Weave Featurizer for molecular graphs - Compatible with scikit-learn models  
* **molgraphconv**: Molecular Graph Convolution Featurizer - Required for GCN model

.. note::
   While scikit-learn models (linear_regression, random_forest_*) can work with any featurizer, the GCN model specifically requires the ``molgraphconv`` featurizer for proper graph-based processing.

Training
--------

The training primitive builds and trains machine learning models on featurized datasets using DeepChem's extensive model library.

.. autofunction:: deepchem_server.core.train.train

Available Models
~~~~~~~~~~~~~~~~

The training primitive supports the following specific model types:

**Scikit-learn Models (wrapped in DeepChem SklearnModel):**

* **linear_regression**: Linear regression for continuous target variables
* **random_forest_classifier**: Random forest for classification tasks  
* **random_forest_regressor**: Random forest for regression tasks

**DeepChem Neural Network Models:**

* **gcn**: Graph Convolutional Network (requires ``molgraphconv`` featurizer)

.. note::
   The GCN model requires PyTorch to be installed and may not be available if torch dependencies are missing. Each model supports different initialization and training parameters - refer to ``deepchem_server.core.model_mappings`` for detailed parameter options.

Inference
---------

The inference primitive runs predictions on new data using previously trained models, handling both featurized and raw input data.

.. autofunction:: deepchem_server.core.inference.infer

Supporting Functions
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: deepchem_server.core.inference._infer_with_featurize

.. autofunction:: deepchem_server.core.inference._infer_without_featurize

Evaluation
----------

The evaluation primitive assesses model performance using various metrics and generates evaluation reports.

.. autofunction:: deepchem_server.core.evaluator.model_evaluator

Supporting Functions
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: deepchem_server.core.evaluator.prc_auc_curve

Molecular Docking
------------------

The docking primitive performs molecular docking between proteins and ligands using AutoDock VINA to predict binding poses and affinities.

.. autofunction:: deepchem_server.core.docking.generate_pose

Available Metrics
~~~~~~~~~~~~~~~~~

The evaluation primitive supports the following metrics:

* **pearson_r2_score**: Pearson correlation coefficient
* **jaccard_score**: Jaccard similarity score
* **prc_auc_score**: Precision-Recall AUC score
* **roc_auc_score**: ROC AUC score
* **rms_score**: Root Mean Square score
* **mae_error**: Mean Absolute Error
* **bedroc_score**: BEDROC score
* **accuracy_score**: Classification accuracy
* **balanced_accuracy_score**: Balanced classification accuracy

Workflow Integration
--------------------

These primitives are designed to work together in typical machine learning workflows:

1. **Data Preparation**: Upload raw data to the datastore
2. **Featurization**: Use ``featurize()`` to transform molecular data into features
3. **Training**: Use ``train()`` to build models on the featurized data
4. **Inference**: Use ``infer()`` to make predictions on new data
5. **Evaluation**: Use ``model_evaluator()`` to assess model performance
6. **Docking**: Use ``generate_pose()`` to predict protein-ligand binding interactions

Example Workflow
~~~~~~~~~~~~~~~~

Here's a typical workflow using all five primitives:

.. code-block:: python

   from deepchem_server.core import feat, train, inference, evaluator, docking
   from deepchem_server.core import config
   from deepchem_server.core.datastore import DiskDataStore
   import tempfile

   # Setup datastore
   datastore = DiskDataStore('profile', 'project', tempfile.mkdtemp())
   config.set_datastore(datastore)

   # 1. Featurize raw data
   dataset_address = feat.featurize(
       dataset_address="raw_data_address",
       featurizer="ecfp",
       output="featurized_dataset",
       dataset_column="smiles",
       label_column="target"
   )

   # 2. Train a model  
   model_address = train.train(
       model_type="random_forest_classifier", 
       dataset_address=dataset_address,
       model_name="my_classification_model"
   )

   # 3. Run inference
   predictions_address = inference.infer(
       model_address=model_address,
       data_address="new_data_address",
       output="predictions.csv",
       dataset_column="smiles"
   )

   # 4. Evaluate the model
   evaluator.model_evaluator(
       dataset_addresses=[dataset_address],
       model_address=model_address,
       metrics=["roc_auc_score", "accuracy_score"],
       output_key="evaluation_results"
   )

   # 5. Perform molecular docking
   docking_address = docking.generate_pose(
       protein_address="protein_address",
       ligand_address="ligand_address",
       output="docking_results"
   ) 