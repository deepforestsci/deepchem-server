# Evaluate

**pyds class:** `Evaluate`
**Server primitive:** `model_evaluator()` in `deepchem_server/core/primitives/evaluator.py`
**Endpoint:** `POST /primitive/evaluate`

## Description

Computes performance metrics for a trained model against one or more featurized datasets. Supports two modes:

1. **Scalar metrics** (`is_metric_plots=False`) — computes one or more numeric scores (ROC-AUC, RMSE, accuracy, etc.) across any number of datasets. Output is a JSON file.
2. **Curve metrics** (`is_metric_plots=True`) — computes a precision-recall curve for a single dataset. Output is a CSV file. Only `prc_auc_curve` is supported in this mode.

## Inputs

| Parameter           | Type        | Required | Default | Description                                                                                                   |
| ------------------- | ----------- | -------- | ------- | ------------------------------------------------------------------------------------------------------------- |
| `dataset_addresses` | `list[str]` | Yes      | —       | **Always a list**, even for a single dataset. Each address must point to a featurized DeepChem `DiskDataset`. |
| `model_address`     | `str`       | Yes      | —       | Datastore address of the trained model.                                                                       |
| `metrics`           | `list[str]` | Yes      | —       | List of metric names to compute. Must match the `is_metric_plots` flag. See valid values below.               |
| `output_key`        | `str`       | Yes      | —       | Short artifact name for the output file. `.json` or `.csv` extension auto-appended based on mode.             |
| `is_metric_plots`   | `bool`      | No       | `False` | Set to `True` to compute curve-based metrics. Enforces single dataset + single metric.                        |

### Valid `metrics` values

**Scalar metrics** (use with `is_metric_plots=False`):

| Value                       | Task type                         |
| --------------------------- | --------------------------------- |
| `"roc_auc_score"`           | Binary classification             |
| `"accuracy_score"`          | Classification                    |
| `"balanced_accuracy_score"` | Classification (class-imbalanced) |
| `"prc_auc_score"`           | Binary classification             |
| `"jaccard_score"`           | Binary classification             |
| `"bedroc_score"`            | Virtual screening / ranking       |
| `"pearson_r2_score"`        | Regression                        |
| `"rms_score"`               | Regression (RMSE)                 |
| `"mae_error"`               | Regression (MAE)                  |

**Curve metrics** (use with `is_metric_plots=True`):

| Value             | Description                                                                                                                                                                   |
| ----------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `"prc_auc_curve"` | Precision-recall curve. Returns `precision`, `recall`, `thresholds` columns. Assumes binary classification with positive-class probabilities in the second prediction column. |

### Constraints

- `is_metric_plots=True` enforces **exactly one metric** and **exactly one dataset**.
- If only plot metrics are passed with `is_metric_plots=False`, raises `ValueError`.
- `prc_auc_curve` uses `y_preds[:, 1]` as the positive-class score — model must output 2-column probabilities.

## Output

**Return key:** `"evaluation_result_address"`
**Return type:** `str`

### Scalar metrics output (JSON)

```json
{
    "deepchem://profile/project/test_set": {
        "roc_auc_score": 0.87,
        "accuracy_score": 0.83
    }
}
```

The top-level keys are the full datastore addresses of the evaluated datasets. Each value is a dict of metric name -> score.

### Curve metrics output (CSV)

| Column       | Description                                                          |
| ------------ | -------------------------------------------------------------------- |
| `precision`  | Precision at each threshold                                          |
| `recall`     | Recall at each threshold                                             |
| `thresholds` | Decision threshold value (last row is `None` per sklearn convention) |
