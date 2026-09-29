# Featurize

**pyds class:** `Featurize`
**Server primitive:** `featurize()` in `deepchem_server/core/primitives/feat.py`
**Endpoint:** `POST /primitive/featurize`

## Description

Converts a raw molecular dataset (CSV or SDF) into a numerical feature representation suitable for machine learning. Applies a molecular featurizer to each molecule, producing a DeepChem `DiskDataset`. Supports both single-core and automatic multi-core featurization for large files (>250 MB), with checkpoint-based restart if interrupted.

## Inputs

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `dataset_address` | `str` | Yes | — | Datastore address of the input dataset. Must be a `.csv` or `.sdf` file. |
| `featurizer` | `str` | Yes | — | Featurization algorithm to apply. See valid values below. |
| `output` | `str` | Yes | — | Short artifact name for the output DiskDataset. Stored as `deepchem://profile/project/<output>`. |
| `dataset_column` | `str` | Yes | — | Column name containing the molecular input (e.g., SMILES strings for CSV). Required for CSV; used by SDF loaders. |
| `feat_kwargs` | `dict` | No | `{}` | Keyword arguments passed to the featurizer constructor (e.g., `{"radius": 2, "size": 1024}` for ECFP). |
| `label_column` | `str` | No | `None` | Column name containing the training labels/targets. **Must be provided for supervised learning datasets** — omitting it produces a dataset with no labels, causing silent failures in training. |
| `n_core` | `int` | No | `os.cpu_count()` | Number of CPU cores for multi-core featurization. |
| `single_core_threshold` | `int` | No | `250` | File size in MB above which multi-core featurization is triggered (when `n_core > 1`). |

### Valid `featurizer` values

| Value | Description |
|-------|-------------|
| `"ecfp"` | Extended Connectivity Fingerprint (Morgan). Fast; good for small molecules. |
| `"graphconv"` | Graph Convolutional features. For GCN models. |
| `"weave"` | Weave/Interaction features. For intermolecular modelling. |
| `"molgraphconv"` | Molecular Graph features. For message-passing networks. |
| `"dummy"` | Identity featurizer. Passes raw input through unchanged. |
| `"grover"` | Pre-trained GROVER transformer features (384-dim). For transfer learning. |
| `"rdkitconformer"` | 3D RDKit conformer features. Requires SDF with 3D coordinates. |
| `"dmpnn"` | Directed Message Passing Neural Network features. |

## Output

**Return key:** `"featurized_file_address"`
**Return type:** `str`

The returned address points to a **DeepChem `DiskDataset` directory** in the datastore containing:
- `X` — feature matrix (numpy array, shape `[n_samples, n_features]`)
- `y` — labels array (empty if `label_column` not provided)
- `ids` — molecule identifiers

A metadata card (`.cdc`) is attached recording the featurizer name and `feat_kwargs` — this is read by `Infer` to auto-featurize raw CSV inputs at inference time.

### Multi-core behavior
- If `file_size > single_core_threshold` AND `n_core > 1`: dataset is split into `n_core` partitions, each featurized in a separate process, then merged.
- Partial results are checkpointed to datastore as `.partial/_checkpoint/part_i_of_n` and cleaned up on success.
- If a run was interrupted, restoring the same `output` key resumes from the last checkpoint.
