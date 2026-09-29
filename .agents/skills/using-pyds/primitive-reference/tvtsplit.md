# TVTSplit (Train-Validation-Test Split)

**pyds class:** `TVTSplit`
**Server primitive:** `train_valid_test_split()` in `deepchem_server/core/primitives/splitter.py`
**Endpoint:** `POST /primitive/train-valid-test-split`

## Description

Splits a featurized dataset into three disjoint subsets — train, validation, and test — according to specified proportions and a splitting strategy. Supports multiple splitting strategies including scaffold-aware splitting for molecular datasets.

## Inputs

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `splitter_type` | `str` | Yes | — | Strategy for splitting. See valid values below. |
| `dataset_address` | `str` | Yes | — | Datastore address of the dataset to split. Supports `DiskDataset`, `NumpyDataset`, and CSV/DataFrame. |
| `frac_train` | `float` | No | `0.8` | Fraction of data assigned to the training set. |
| `frac_valid` | `float` | No | `0.1` | Fraction of data assigned to the validation set. |
| `frac_test` | `float` | No | `0.1` | Fraction of data assigned to the test set. |

Fractions should sum to 1.0. No validation is performed — if they don't sum to 1.0, data will be silently lost or duplicated.

### Valid `splitter_type` values

| Value | Description |
|-------|-------------|
| `"random"` | Randomly shuffles all samples before splitting. Default choice for most tasks. |
| `"scaffold"` | Groups molecules by Murcko scaffold, then splits — ensures no scaffold appears in both train and test. Recommended for unbiased molecular property models. |
| `"index"` | Splits by index order (no shuffle). Useful when samples have a meaningful ordering. |
| `"random_stratified"` | Stratified random split that preserves class-label proportions across splits. Use for imbalanced classification tasks. |

## Output

**Return key:** `"train_valid_test_split_results_address"`
**Return type:** `list[str]` — exactly 3 elements

```python
[train_address, valid_address, test_address]  # indices 0, 1, 2
```

Each address points to a dataset of the same type as the input (DiskDataset or CSV), stored as:
- `<key>_train`
- `<key>_valid`
- `<key>_test`

The featurizer metadata from the parent dataset is preserved in each split's data card.

### Notes
- `splitter_type` only applies scaffold/stratified logic for `DiskDataset`/`NumpyDataset`. **CSV datasets always use random shuffling** regardless of `splitter_type`.
- The return value is a list — always index with `[0]` (train), `[1]` (valid), `[2]` (test).

---

# k_fold_split

**pyds class:** Not yet exposed via pyds client (server-side only)
**Server primitive:** `k_fold_split()` in `deepchem_server/core/primitives/splitter.py`

## Description

Generates k train/test fold pairs for cross-validation. Each fold uses approximately 1/k of the data as the test set and the remaining as training.

## Inputs

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `splitter_type` | `str` | Yes | Same valid values as TVTSplit. |
| `dataset_address` | `str` | Yes | Dataset to fold. |
| `k` | `int` | Yes | Number of folds. |

## Output

`list[tuple[str, str]]` — k tuples of `(train_address, test_address)`.

Each fold is stored as `<key>_train{i}` and `<key>_test{i}` where `i` is 0-indexed.
