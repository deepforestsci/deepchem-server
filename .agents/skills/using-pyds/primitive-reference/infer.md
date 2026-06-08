# Infer

**pyds class:** `Infer`
**Server primitive:** `infer()` in `deepchem_server/core/primitives/inference.py`
**Endpoint:** `POST /primitive/infer`

## Description

Runs a trained model on new data and writes predictions to a CSV file. Supports two data paths:

1. **Auto-featurization path** — input is a raw CSV; the server looks up the featurizer from the model's `ModelCard -> training dataset card` chain and featurizes on-the-fly in chunks.
2. **Direct prediction path** — input is an already-featurized `DiskDataset`; predictions run directly without re-featurizing.

Data is processed in streaming chunks (`shard_size` rows at a time) to avoid loading large datasets into memory.

## Inputs

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `model_address` | `str` | Yes | — | Datastore address of the trained model (output of `Train`). |
| `data_address` | `str` | Yes | — | Datastore address of the inference dataset. **Must be a datastore address** — not a local file path. Upload raw CSV first with `Data.upload_data()`. |
| `output` | `str` | Yes | — | Short artifact name for the output predictions CSV. `.csv` extension auto-appended if missing. |
| `dataset_column` | `str` | No | `None` | Column name containing the raw molecular inputs (e.g., SMILES). **Required when `data_address` points to a raw CSV** — triggers auto-featurization. Ignored if data is already featurized. |
| `shard_size` | `int` | No | `8192` | Number of rows processed per chunk. Lower values reduce memory use; higher values improve throughput. |
| `threshold` | `float` | No | `None` | Decision threshold for binary classification. If provided, appends a `binarized_preds` column: rows where the predicted probability exceeds the threshold are marked `1`, others `0`. Applied to the last prediction column only. |

### Auto-featurization trigger condition
The server uses auto-featurization when **both** conditions are true:
- The data card has no `featurizer` field (raw data, not a DiskDataset)
- `data_address` ends with `.csv`

If these conditions are met, `dataset_column` becomes mandatory.

## Output

**Return key:** `"inference_results_address"`
**Return type:** `str`

The returned address points to a **CSV file** in the datastore.

### Output CSV structure

| Column | Description |
|--------|-------------|
| `X` | Raw inputs (SMILES strings for auto-featurized data) or sample IDs (for pre-featurized data) |
| `y_preds` | Predictions (single-task models) |
| `y1_preds`, `y2_preds`, … | Predictions per task (multi-task models) |
| `binarized_preds` | Binary predictions (only present if `threshold` is set) |

### Notes
- For multi-task models, each task gets its own column: `y1_preds`, `y2_preds`, etc.
- Threshold binarization applies to the last prediction column only — appropriate for binary classification where the last column is the positive-class probability.
- The featurizer lookup chain (Model card → training dataset card) must remain intact for auto-featurization to work. If the training dataset was deleted, pass a pre-featurized dataset instead.
