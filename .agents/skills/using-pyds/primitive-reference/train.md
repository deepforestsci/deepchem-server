# Train

**pyds class:** `Train`
**Server primitive:** `train()` in `deepchem_server/core/primitives/train.py`
**Endpoint:** `POST /primitive/train`

## Description

Trains a machine learning model on a featurized DeepChem dataset. Accepts both sklearn-style models (random forest, linear regression) and deep learning models (GCN via PyTorch). Stores the trained model and a metadata card (`ModelCard`) in the datastore — the card records the training dataset address so that `Infer` can automatically look up the correct featurizer later.

## Inputs

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `dataset_address` | `str` | Yes | — | Datastore address of the featurized training dataset. Must be a DeepChem `DiskDataset` (i.e., output of `Featurize`). |
| `model_type` | `str` | Yes | — | Model architecture identifier. See valid values below. |
| `model_name` | `str` | Yes | — | Short artifact name for the saved model. Stored as `deepchem://profile/project/<model_name>`. Raises `FileExistsError` if already exists. |
| `init_kwargs` | `dict` | No | `{}` | Keyword arguments for the model constructor (e.g., `{"n_estimators": 100}` for random forest, `{"graph_conv_layers": [64, 64]}` for GCN). |
| `train_kwargs` | `dict` | No | `{}` | Keyword arguments for `model.fit()`. For deep learning: `{"nb_epoch": 50}`. |

### Valid `model_type` values

| Value | Type | Notes |
|-------|------|-------|
| `"random_forest_classifier"` | sklearn | Binary or multi-class classification. |
| `"random_forest_regressor"` | sklearn | Continuous output regression. |
| `"linear_regression"` | sklearn | Linear regression. |
| `"gcn"` | PyTorch (DeepChem) | Graph Convolutional Network. Requires `n_tasks` in `init_kwargs`. Progress logged during training. |

## Output

**Return key:** `"trained_model_address"`
**Return type:** `str`

The returned address points to the **serialized model directory** in the datastore.

A `ModelCard` metadata file (`.cmc`) is stored alongside the model containing:
- `model_type` — the model type string
- `train_dataset_address` — address of the training dataset (**critical**: used by `Infer` to look up the featurizer for auto-featurization of raw CSV inputs)
- `init_kwargs`, `train_kwargs` — parameters used, for reproducibility

### Important notes
- The `train_dataset_address` stored in `ModelCard` is the featurizer lookup chain for inference. If the training dataset or its card is deleted, inference on raw CSV will fail.
- For GCN, `nb_epoch` in `train_kwargs` controls training length (default: 10 if not specified). For production use, pass higher values.
