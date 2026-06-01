---
name: using-pyds
description: Use when writing Python code that interacts with the DeepChem Server via the pyds client library — uploading/downloading data, running primitives (featurize, train, infer, evaluate, split, partition, docking, DEL denoising, PDB cleaning, RBFE), or chaining primitive outputs as inputs to build computational chemistry workflows.
---

# pyds API Reference

## Python Environment

<!-- PYDS_PYTHON_BOOTSTRAP_START -->

**First use only — Python environment:** `PYDS_PYTHON` is not yet set for this project. Ask the user for the default Python interpreter path (e.g. `conda run -n deepchem-server-env python` or an absolute path). Once confirmed, edit this skill file to:

1. Remove the `PYDS_PYTHON_BOOTSTRAP_START` … `PYDS_PYTHON_BOOTSTRAP_END` block entirely
2. Add a single line immediately after the `## Python Environment` heading: `PYDS_PYTHON: conda run -n deepchem-server-env python`
 <!-- PYDS_PYTHON_BOOTSTRAP_END -->

## Session Setup

Do this at the start of every pyds script or session, before any other operations.

### Step 1 — Imports

```python
from pyds import Settings, Data, Featurize, TVTSplit, Train, Evaluate, Infer
from pyds import Docking, DelDenoise, Partition, BaseClient
from pyds.primitives import PdbClean  # NOT in top-level pyds; must import from pyds.primitives
```

### Step 2 — Load or bootstrap settings

Always initialize with no arguments. `Settings()` loads `.pyds.settings.json` from the **current working directory** automatically — creates an empty file if missing, treats parse errors as unconfigured (no exception raised).

```python
settings = Settings()

if not settings.is_configured():
    # Ask the user (in conversation) for the required values, then apply.
    # Each set_* method saves to .pyds.settings.json immediately,
    # so subsequent sessions load without prompting.
    settings.set_profile("<ask user: profile name>")
    settings.set_project("<ask user: project name>")
    # base_url defaults to http://localhost:8000 — ask only if the user needs a different URL
    settings.set_base_url("<ask user: server URL, or skip to keep http://localhost:8000>")
```

`is_configured()` returns `True` only when both `profile` and `project` are non-null. A missing file, empty file, or corrupt JSON all result in `is_configured() == False` and trigger the prompt above.

### Step 3 — Health check

Run once after settings are confirmed. Uses no profile/project — works with any initialized settings.

```python
health = BaseClient(settings=settings).healthcheck()
# Returns {"status": ..., "version": ...}
# Raises an exception if the server is unreachable — stop and report the error; do not proceed.
print(health)
```

Pass `settings` to every subsequent client: `Data(settings=settings)`, `Train(settings=settings)`, etc.

## Address System

Every datastore object is identified by a **DeepChem address**. Primitive outputs return addresses; those addresses become inputs to the next primitive.

**Canonical format:** `deepchem://profile/project/key`

- `profile` — first path component (matches `Settings.profile`)
- `project` — second path component (matches `Settings.project`)
- `key` — everything after `profile/project/` — may contain `/` (e.g. `models/rf_v1`, `data/featurized/ecfp`)

**Short form** (prefix omitted) is also accepted as input to any primitive: `profile/project/key`. The server normalizes it to the canonical form internally. Minimum 3 path components are always required — a bare filename like `"zinc.csv"` is invalid.

**Returned addresses always include the `deepchem://` prefix.** Pass them directly to the next primitive without modification.

```
deepchem://alice/qsar_project/data/molecules.csv      ← returned by upload_data
deepchem://alice/qsar_project/features/ecfp           ← returned by Featurize
deepchem://alice/qsar_project/models/rf_v1            ← returned by Train
```

Each primitive returns its result under a unique key — see the Primitives table below.

## Data Operations

```python
client = Data(settings=settings)

# Upload a local file
result = client.upload_data(file_path="molecules.csv", description="optional")
addr = result["dataset_address"]
# addr is a full deepchem:// address, e.g. "deepchem://alice/qsar/data/molecules.csv"

# List all files in the project
files = client.list_data()  # returns list of deepchem:// address strings

# Download to local path — pass the full address as returned by upload or a primitive
local_path = client.download_data(address=addr, destination_path="/tmp/out.csv")
```

All methods accept optional `profile_name`, `project_name`, `backend` overrides.

## Primitives — Signatures and Return Keys

All primitives: `client = ClassName(settings=settings)` then `result = client.run(...)`.

| Primitive  | Class        | Required params                                                                                               | Key optional params                                                                        | Return key                               | Return type                          |
| ---------- | ------------ | ------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ | ---------------------------------------- | ------------------------------------ |
| Featurize  | `Featurize`  | `dataset_address, featurizer, output, dataset_column`                                                         | `feat_kwargs={}, label_column`                                                             | `featurized_file_address`                | `str`                                |
| Train      | `Train`      | `dataset_address, model_type, model_name`                                                                     | `init_kwargs={}, train_kwargs={}`                                                          | `trained_model_address`                  | `str`                                |
| Infer      | `Infer`      | `model_address, data_address, output`                                                                         | `dataset_column, shard_size=8192, threshold`                                               | `inference_results_address`              | `str`                                |
| Evaluate   | `Evaluate`   | `dataset_addresses` (list), `model_address, metrics` (list), `output_key`                                     | `is_metric_plots=False`                                                                    | `evaluation_result_address`              | `str`                                |
| TVTSplit   | `TVTSplit`   | `splitter_type, dataset_address`                                                                              | `frac_train=0.8, frac_valid=0.1, frac_test=0.1`                                            | `train_valid_test_split_results_address` | `list[str]` → `[train, valid, test]` |
| Partition  | `Partition`  | `dataset_address`                                                                                             | `n_partition=4, shuffle=False`                                                             | `partitioned_dataset_addresses`          | `list[str]`                          |
| Docking    | `Docking`    | `protein_address, ligand_address, output`                                                                     | `exhaustiveness=10, num_modes=9`                                                           | `docking_results_address`                | `str`                                |
| DelDenoise | `DelDenoise` | `dataset_address, output_key`                                                                                 | `strategy="unified", control_cols, target_cols, add_hit_labels=False, hit_percentile=90.0` | `denoised_dataset_address`               | `str`                                |
| PdbClean   | `PdbClean`   | `gene_name` (human gene, e.g. `"MAP2K1"`), `entry_id` (UniProt accession, e.g. `"P02144"`), `scan_id, output` | `min_res_val=2.5`                                                                          | `pdb_clean_results_address`              | `str`                                |

## Supported Values

**`featurizer`:** `ecfp`, `graphconv`, `weave`, `molgraphconv`, `dummy`, `grover`, `rdkitconformer`, `dmpnn`

**`model_type`:** `linear_regression`, `random_forest_regressor`, `random_forest_classifier`, `gcn`

**`splitter_type`:** `random`, `index`, `scaffold`, `random_stratified`

**`metrics`:** `pearson_r2_score`, `roc_auc_score`, `accuracy_score`, `balanced_accuracy_score`, `rms_score`, `mae_error`, `prc_auc_score`, `jaccard_score`, `bedroc_score`, `prc_auc_curve`

## Data Type Flow

```
Raw CSV/SDF  →  Featurize  →  DiskDataset (directory)
DiskDataset  →  TVTSplit   →  [train_addr, valid_addr, test_addr]
DiskDataset  →  Train      →  model directory
model + data →  Infer      →  predictions CSV
model + data →  Evaluate   →  metrics JSON (or CSV for prc_auc_curve)
PDB file     →  PdbClean   →  cleaned PDB  →  Docking (+ SDF ligand)  →  results JSON
DEL CSV      →  DelDenoise →  scored CSV
```

**Note:** `Infer.data_address` must be a datastore address (not a local file path). Upload raw CSV first via `Data.upload_data()`, then pass `dataset_address`.

## `output` Parameter

All primitives that accept `output` take a **short artifact name** (not a full path or address). The server stores the result under `profile/project/<type>/<output>`. Examples: `"ecfp_features"`, `"rf_model_v1"`, `"dock_result"`.

## Primitive Deep Reference

Detailed documentation for each primitive (inputs, output format, algorithm notes, edge cases) lives in `primitive-reference/`:

| File                                                                       | Primitive                                                 |
| -------------------------------------------------------------------------- | --------------------------------------------------------- |
| [`primitive-reference/featurize.md`](primitive-reference/featurize.md)     | `Featurize` — CSV/SDF → DiskDataset                       |
| [`primitive-reference/train.md`](primitive-reference/train.md)             | `Train` — DiskDataset → Model                             |
| [`primitive-reference/infer.md`](primitive-reference/infer.md)             | `Infer` — Model + Data → predictions CSV                  |
| [`primitive-reference/evaluate.md`](primitive-reference/evaluate.md)       | `Evaluate` — Model + Datasets → metrics JSON/CSV          |
| [`primitive-reference/tvtsplit.md`](primitive-reference/tvtsplit.md)       | `TVTSplit` — Dataset → [train, valid, test]               |
| [`primitive-reference/partition.md`](primitive-reference/partition.md)     | `Partition` — Dataset → N partitions                      |
| [`primitive-reference/docking.md`](primitive-reference/docking.md)         | `Docking` — Protein + Ligand → poses + scores             |
| [`primitive-reference/del_denoise.md`](primitive-reference/del_denoise.md) | `DelDenoise` — DEL counts CSV → enrichment scores         |
| [`primitive-reference/pdb_clean.md`](primitive-reference/pdb_clean.md)     | `PdbClean` — UniProt ID → cleaned PDB(s) + summary JSON   |
| [`primitive-reference/rbfe.md`](primitive-reference/rbfe.md)               | `run_rbfe` / `collate_rbfe_results` — FEP ΔΔG simulations |

Read the relevant file when you need parameter-level detail, output field names, or non-obvious constraints for a specific primitive.

## Common Mistakes

- **`PdbClean` import:** Use `from pyds.primitives import PdbClean` — it is not in `from pyds import *`.
- **`Evaluate.dataset_addresses` is always a list:** `[test_addr]` not `test_addr`.
- **`TVTSplit` return value is a list:** `result["train_valid_test_split_results_address"]` → index `[0]` train, `[1]` valid, `[2]` test.
- **Featurizing supervised data:** Always pass `label_column` when the dataset has a target column; omitting it produces a dataset with no labels.
- **`Infer` on raw CSV:** Pass `dataset_column` (the SMILES column name) so the server auto-featurizes using the model's original featurizer.
- **Return key spelling:** Each primitive has a unique key name — never assume `"result_address"`. Use the table above.
- **Bare filenames as addresses are invalid:** `"zinc.csv"` will raise `ValueError: Invalid deepchem address`. Always use the full `deepchem://profile/project/key` form (or at minimum `profile/project/key`). Use the address returned by the previous primitive — don't construct addresses by hand.
