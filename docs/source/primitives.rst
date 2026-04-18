Overview
--------

The core primitives that form the backbone of DeepChem Server's machine learning workflow.
These primitives provide the essential functionality for molecular machine learning pipelines.
Currently, Deepchem Server provides the following primitives: (Other primitives are planned to be added soon)

* **Featurize**: Transform raw molecular data into machine learning features
* **Partition**: Split a dataset into multiple datastore-backed partitions
* **Train Valid Test Split**: Split the dataset into training, validation, and test sets
* **Train**: Build and train machine learning models on featurized datasets
* **Inference**: Run predictions on new data using trained models
* **Evaluation**: Assess model performance using various metrics
* **Docking**: Perform molecular docking to predict protein-ligand binding poses
* **DEL Denoise**: Score DEL screening data to identify strong binders
* **PDB clean**: Prepare protein structures (PDBFixer/OpenMM: heterogens, water, chains, hydrogens)
* **Ligand prep**: Build 3D ligand conformers from SMILES (RDKit) and store as SDF in the datastore
* **Filter promiscuous targets**: Identify and remove promiscuous binding targets from per-ligand docking scan results

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

Partition
---------

The partition primitive splits a dataset into multiple smaller datasets and uploads each partition back to the datastore. Use it when you need to:

* **Distribute workloads** — break a large dataset into chunks so featurization or inference can run in parallel.
* **Reduce memory pressure** — process data in manageable pieces when the full dataset does not fit in memory.
* **Stage training pipelines** — create fixed data splits before training begins.

**Supported dataset types:**

* DeepChem ``DiskDataset`` — supports optional shuffling before partitioning.
* CSV ``DataFrame`` — partitions rows sequentially (shuffling is not supported).

After partitioning, the parent dataset's datacard is updated with the number of partitions (n_partition) so downstream primitives can discover the splits.

.. autofunction:: deepchem_server.core.primitives.partition.partition

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

**Key Features:**

* Generates protein-ligand binding poses using AutoDock VINA
* Supports both PDB and PDBQT output formats
* Automatically splits PDBQT files for multiple binding modes
* Returns DeepChem addresses to all generated files


.. autofunction:: deepchem_server.core.docking.generate_pose

Supporting Functions
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: deepchem_server.core.docking.split_pdbqt_docked_ligands
   :no-index:

PDB cleaning
------------

The PDB cleaning primitive loads a structure from the datastore, optionally removes chains, 
heterogens, and water, adds missing hydrogens at a chosen pH using PDBFixer and OpenMM, 
and writes a cleaned PDB back to the datastore.

**Dependencies:** PDBFixer and OpenMM must be installed in the execution environment.

.. autofunction:: deepchem_server.core.primitives.pdb_clean.pdb_clean

Ligand preparation
------------------

The ligand preparation primitive converts a SMILES string to a 3D molecule with RDKit (ETKDG),
serializes it as SDF, and uploads the file to the datastore.

**Dependencies:** RDKit must be installed in the execution environment.

.. autofunction:: deepchem_server.core.primitives.ligand_prep.ligand_prep

Filter promiscuous targets
--------------------------

The filter promiscuous targets primitive identifies and removes promiscuous binding targets from
per-ligand docking scan results. A promiscuous target is a gene that appears in the top-M% of
results across N or more ligands, suggesting non-selective binding.

**Input Requirements:**

* Per-ligand scan result CSV files from docking workflows
* Each CSV must contain a `gene_name` column
* Rows should be sorted from best (top) to worst docking score

**Threshold Parameters:**

The primitive accepts a list of [m, n] threshold pairs:

* **m**: Percentile cutoff (0-100) — top m% of results to consider
* **n**: Minimum occurrence count — gene must appear in at least n ligands' top results

**Outputs:**

* JSON file mapping thresholds to promiscuous gene lists
* Filtered CSV files with promiscuous targets removed

.. autofunction:: deepchem_server.core.primitives.filter_promiscuous_targets.filter_promiscuous_targets

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

DEL Denoise
-----------

The DEL Denoise primitive scores DEL screening data to identify compounds that 
are strongly enriched in the target selection relative to background noise.

**Key Features:**

* Supports two scoring strategies: Poisson-based (``unified``) and z-score-based (``non_unified``)
* Optional collapsing of trisynthon rows into pairwise disynthon combinations

**Scoring Strategies:**

* **unified**: Computes a Poisson confidence-interval enrichment ratio
  (target lower bound / control upper bound) across all replicates simultaneously.
* **non_unified**: Sums replicates then computes a z-score for each compound
  independently in the target and control.

.. autofunction:: deepchem_server.core.primitives.del_denoising.del_denoise

Supporting Functions
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: deepchem_server.core.primitives.del_denoising.poissfit
   :no-index:

.. autofunction:: deepchem_server.core.primitives.del_denoising.get_enrichment_ratio
   :no-index:

.. autofunction:: deepchem_server.core.primitives.del_denoising.calculate_poisson_enrichment
   :no-index:

.. autofunction:: deepchem_server.core.primitives.del_denoising.calculate_normalized_enrichment_score
   :no-index:

.. autofunction:: deepchem_server.core.primitives.del_denoising.calculate_hit_threshold
   :no-index:

.. autofunction:: deepchem_server.core.primitives.del_denoising.collapse_to_disynthons
   :no-index:

Workflow Integration
--------------------

These primitives are designed to work together in typical machine learning workflows:

1. **Data Preparation**: Upload raw data to the datastore
2. **Partition (Optional)**: Use ``partition()`` to split large datasets for parallel or staged workflows
3. **Featurization**: Use ``featurize()`` to transform molecular data into features
4. **Training**: Use ``train()`` to build models on the featurized data
5. **Inference**: Use ``infer()`` to make predictions on new data
6. **Evaluation**: Use ``model_evaluator()`` to assess model performance
7. **Docking**: Use ``generate_pose()`` to predict protein-ligand binding interactions
8. **PDB clean** (optional): Use ``pdb_clean()`` to standardize a receptor PDB before docking or simulation
9. **Ligand prep** (optional): Use ``ligand_prep()`` to generate a 3D SDF for a ligand from SMILES
10. **Filter promiscuous targets** (optional): Use ``filter_promiscuous_targets()`` to remove non-selective targets from per-ligand docking results

Example Workflow
~~~~~~~~~~~~~~~~

Here's a typical workflow using all five primitives:

.. code-block:: python

   from deepchem_server.core import feat, train, inference, evaluator, docking
   from deepchem_server.core.primitives.pdb_clean import pdb_clean
   from deepchem_server.core.primitives.ligand_prep import ligand_prep
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
       protein_address=cleaned_protein,
       ligand_address=ligand_sdf,
       output="docking_results"
   ) 