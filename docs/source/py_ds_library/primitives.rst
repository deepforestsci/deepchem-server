Primitives
==========

The primitives system provides a modular approach to computation tasks. All primitives inherit from the base Primitive class and implement a standardized interface.

Primitive Base Class
--------------------

The Primitive class is an abstract base class that all computation primitives inherit from.

.. autoclass:: pyds.primitives.base.Primitive
   :members:
   :undoc-members:

   .. note::
      This is an abstract base class. Concrete implementations must override the ``run()`` method.

Featurize Primitive
-------------------

The Featurize primitive transforms raw molecular data into machine learning features using DeepChem's featurizers.

.. autoclass:: pyds.primitives.featurize.Featurize
   :members:
   :undoc-members:

   **Common Featurizers:**

   * **ECFP**: Extended Connectivity Fingerprints (Circular Fingerprints)
   * **GraphConv**: Graph Convolution Featurizer
   * **Weave**: Weave Featurizer for molecular graphs
   * **MolGraphConv**: Molecular Graph Convolution Featurizer (required for GCN models)

   **Example Usage:**

   .. code-block:: python

      featurize = Featurize(settings)
      response = featurize.run(
          dataset_address="dataset_address",
          featurizer="ECFP",
          output="featurized_dataset",
          dataset_column="smiles",
          feat_kwargs={"radius": 2, "size": 1024},
          label_column="target"
      )

Train Primitive
---------------

The Train primitive builds and trains machine learning models on featurized datasets.

.. autoclass:: pyds.primitives.train.Train
   :members:
   :undoc-members:

   **Supported Model Types:**

   * **Scikit-learn Models**: ``linear_regression``, ``random_forest_classifier``, ``random_forest_regressor``
   * **DeepChem Neural Networks**: ``gcn`` (Graph Convolutional Network)

   **Example Usage:**

   .. code-block:: python

      train = Train(settings)
      response = train.run(
          dataset_address="featurized_dataset_address",
          model_type="random_forest_classifier",
          model_name="my_classifier",
          init_kwargs={"n_estimators": 100},
          train_kwargs={"nb_epoch": 10}
      )

Evaluate Primitive
------------------

The Evaluate primitive assesses model performance using various metrics and generates evaluation reports.

.. autoclass:: pyds.primitives.evaluate.Evaluate
   :members:
   :undoc-members:

   **Available Metrics:**

   * **Classification**: ``roc_auc_score``, ``accuracy_score``, ``balanced_accuracy_score``, ``prc_auc_score``
   * **Regression**: ``pearson_r2_score``, ``rms_score``, ``mae_error``
   * **Other**: ``jaccard_score``, ``bedroc_score``

   **Example Usage:**

   .. code-block:: python

      evaluate = Evaluate(settings)
      response = evaluate.run(
          dataset_addresses=["test_dataset_address"],
          model_address="trained_model_address",
          metrics=["roc_auc_score", "accuracy_score"],
          output_key="evaluation_results",
          is_metric_plots=True
      )

Infer Primitive
---------------

The Infer primitive runs predictions on new data using previously trained models.

.. autoclass:: pyds.primitives.infer.Infer
   :members:
   :undoc-members:

   **Key Features:**

   * Handles both featurized and raw input data
   * Supports classification thresholding
   * Configurable shard size for large datasets

   **Example Usage:**

   .. code-block:: python

      infer = Infer(settings)
      response = infer.run(
          model_address="trained_model_address",
          data_address="new_data_address",
          output="predictions",
          dataset_column="smiles",
          threshold=0.5
      )

TVTSplit Primitive
------------------

The TVTSplit primitive performs train-validation-test splitting on datasets using various splitter types.

.. autoclass:: pyds.primitives.splitter.TVTSplit
   :members:
   :undoc-members:

   **Common Splitter Types:**

   * **RandomSplitter**: Random data splitting
   * **ScaffoldSplitter**: Split based on molecular scaffolds
   * **ButinaSplitter**: Cluster-based splitting using molecular fingerprints

   **Example Usage:**

   .. code-block:: python

      split = TVTSplit(settings)
      response = split.run(
          splitter_type="RandomSplitter",
          dataset_address="dataset_address",
          frac_train=0.8,
          frac_valid=0.1,
          frac_test=0.1
      )
