# DEL Denoising and ML Pipelines

DNA-Encoded Library (DEL) screening identifies candidate binders by coupling each compound to a DNA barcode, incubating the library against a protein target, and sequencing the surviving barcodes. The resulting raw sequencing counts are corrupted by PCR amplification bias and stochastic noise, so they cannot be used directly as ML labels. The `del_denoise` primitive converts raw counts into enrichment scores — statistically principled labels that downstream models can learn from.

## Choose Your Workflow

Two orthogonal choices determine which pipeline to run:

| | **Unified** (Poisson CI) | **Non-Unified** (z-score) |
|---|---|---|
| **Trisynthon** | [Workflow A](#workflow-a-trisynthon--unified) — regression or classification; best overall accuracy | [Workflow B](#workflow-b-trisynthon--non-unified) — two-model classification; best for classification tasks |
| **Disynthon** | ⚠️ [Not recommended](#warning-disynthon--unified) | [Workflow C](#workflow-c-disynthon--non-unified) — same as B; 2–4× faster training on large libraries |

**Trisynthon** keeps each compound as its three-fragment SMILES. **Disynthon** collapses all three-fragment compounds into their three pairwise two-fragment combinations, aggregating counts — amplifying signal at the cost of structural information.

**Unified** scores target and control replicates jointly via Poisson confidence intervals, producing one enrichment value per compound. **Non-unified** scores them independently via z-score normalization, producing separate target and control labels that feed two models trained in parallel.

---

## Prerequisites

### Input CSV format

The raw CSV must have these columns. Additional columns are preserved but ignored.

```
smiles,smiles_a,smiles_b,smiles_c,seq_target_1,seq_target_2,seq_target_3,seq_matrix_1,seq_matrix_2,seq_matrix_3
```

| Column | Description |
|--------|-------------|
| `smiles` | Full trisynthon SMILES for the compound |
| `smiles_a`, `smiles_b`, `smiles_c` | Fragment SMILES for synthons A, B, C — required only for disynthon collapse |
| `seq_target_1/2/3` | Sequencing read counts from target selection replicates |
| `seq_matrix_1/2/3` | Sequencing read counts from control (matrix) selection replicates |

### Session setup

```python
from pyds import Settings, Data, Featurize, TVTSplit, Train, Evaluate, Infer, DelDenoise
from pyds.base import BaseClient

settings = Settings()
if not settings.is_configured():
    settings.set_profile("your_profile")
    settings.set_project("your_project")
    # base_url defaults to http://localhost:8000
    settings.set_base_url("http://localhost:8000")

BaseClient(settings=settings).healthcheck()  # raises if server unreachable

data_client   = Data(settings=settings)
denoiser      = DelDenoise(settings=settings)
featurizer    = Featurize(settings=settings)
splitter      = TVTSplit(settings=settings)
trainer       = Train(settings=settings)
evaluator     = Evaluate(settings=settings)
infer_client  = Infer(settings=settings)
```

---

## Workflow A: Trisynthon + Unified

**Output:** `Poisson_Enrichment` — one float per compound representing `target_lower_CI / control_upper_CI`. Values above 1.0 indicate enrichment above background.

**Use for:** regression (predict enrichment directly) or classification (threshold at a percentile with `add_hit_labels=True`).

```python
# 1. Upload raw CSV
upload   = data_client.upload_data(file_path="del_screen.csv")
raw_addr = upload["dataset_address"]

# 2. Denoise — unified Poisson scoring
denoise_result = denoiser.run(
    dataset_address=raw_addr,
    output_key="denoised_unified",
    strategy="unified",
    # add_hit_labels=True, hit_percentile=90.0  ← uncomment for classification
)
denoised_addr = denoise_result["denoised_dataset_address"]

# 3. Featurize
# For regression: label_column="Poisson_Enrichment"
# For classification: label_column="hits"  (requires add_hit_labels=True above)
feat_result = featurizer.run(
    dataset_address=denoised_addr,
    featurizer="ecfp",
    output="del_ecfp",
    dataset_column="smiles",
    label_column="Poisson_Enrichment",
)
feat_addr = feat_result["featurized_file_address"]

# 4. Split
split_result = splitter.run(
    dataset_address=feat_addr,
    splitter_type="random",
    frac_train=0.8,
    frac_valid=0.1,
    frac_test=0.1,
)
train_addr, valid_addr, test_addr = split_result["train_valid_test_split_results_address"]

# 5. Train
# For regression: model_type="random_forest_regressor"
# For classification: model_type="random_forest_classifier"
train_result = trainer.run(
    dataset_address=train_addr,
    model_type="random_forest_regressor",
    model_name="del_rf_unified",
)
model_addr = train_result["trained_model_address"]

# 6. Evaluate
evaluator.run(
    dataset_addresses=[train_addr, valid_addr, test_addr],
    model_address=model_addr,
    metrics=["rms_score"],   # regression; use ["roc_auc_score"] for classification
    output_key="del_eval_unified",
)
```

---

## Workflow B: Trisynthon + Non-Unified

**Output:** `target_hits` and `control_hits` — binary (0/1) hit labels derived independently from target z-scores and control z-scores at the 90th percentile.

**Use for:** classification tasks. Two models are trained — one on the target, one on the control — and their predictions are combined at inference time.

```python
# 1. Upload raw CSV
upload   = data_client.upload_data(file_path="del_screen.csv")
raw_addr = upload["dataset_address"]

# 2. Denoise — non-unified z-score, with hit labels
denoise_result = denoiser.run(
    dataset_address=raw_addr,
    output_key="denoised_non_unified",
    strategy="non_unified",
    add_hit_labels=True,
    hit_percentile=90.0,   # top 10% are hits
)
denoised_addr = denoise_result["denoised_dataset_address"]

# 3–5. Target model: featurize on target_hits, split, train, infer
feat_target = featurizer.run(
    dataset_address=denoised_addr,
    featurizer="ecfp",
    output="del_ecfp_target",
    dataset_column="smiles",
    label_column="target_hits",
)
feat_target_addr = feat_target["featurized_file_address"]

split_target = splitter.run(
    dataset_address=feat_target_addr,
    splitter_type="random",
    frac_train=0.8, frac_valid=0.1, frac_test=0.1,
)
train_t, valid_t, test_t = split_target["train_valid_test_split_results_address"]

target_model = trainer.run(
    dataset_address=train_t,
    model_type="random_forest_classifier",
    model_name="del_target_model",
)
target_model_addr = target_model["trained_model_address"]

target_infer = infer_client.run(
    model_address=target_model_addr,
    data_address=test_t,
    output="del_target_preds",
)
target_infer_addr = target_infer["inference_results_address"]

# 3–5. Control model: same steps on control_hits
feat_control = featurizer.run(
    dataset_address=denoised_addr,
    featurizer="ecfp",
    output="del_ecfp_control",
    dataset_column="smiles",
    label_column="control_hits",
)
feat_control_addr = feat_control["featurized_file_address"]

split_control = splitter.run(
    dataset_address=feat_control_addr,
    splitter_type="random",
    frac_train=0.8, frac_valid=0.1, frac_test=0.1,
)
train_c, valid_c, test_c = split_control["train_valid_test_split_results_address"]

control_model = trainer.run(
    dataset_address=train_c,
    model_type="random_forest_classifier",
    model_name="del_control_model",
)
control_model_addr = control_model["trained_model_address"]

control_infer = infer_client.run(
    model_address=control_model_addr,
    data_address=test_c,
    output="del_control_preds",
)
control_infer_addr = control_infer["inference_results_address"]

# 6. Combine predictions: Target AND NOT Control
# p_combined = p_target × (1 − p_control)
# y2_preds is the positive-class probability column from the classifier
import pandas as pd

data_client.download_data(address=target_infer_addr, destination_path="/tmp/target_preds.csv")
data_client.download_data(address=control_infer_addr, destination_path="/tmp/control_preds.csv")

target_df  = pd.read_csv("/tmp/target_preds.csv")[["X", "y2_preds"]].rename(columns={"y2_preds": "p_target"})
control_df = pd.read_csv("/tmp/control_preds.csv")[["X", "y2_preds"]].rename(columns={"y2_preds": "p_control"})

combined = target_df.merge(control_df, on="X")
combined["p_combined"] = combined["p_target"] * (1 - combined["p_control"])
# Rank compounds by p_combined descending; higher = stronger predicted binder
combined.sort_values("p_combined", ascending=False).to_csv("del_ranked_hits.csv", index=False)
```

---

## Workflow C: Disynthon + Non-Unified

Identical to Workflow B, but `del_denoise` first collapses each trisynthon into its three pairwise disynthon combinations (AB, AC, BC), summing counts across all trisynthons that share a pair. The output SMILES column is `disynthons` (two-fragment canonical SMILES joined by `.`).

**Use for:** very large DEL assays where trisynthon-scale training is computationally prohibitive. Expect 2–4× faster featurization and training at the cost of slightly lower model accuracy vs Workflow B.

```python
# 1. Upload raw CSV (must have smiles_a, smiles_b, smiles_c columns)
upload   = data_client.upload_data(file_path="del_screen.csv")
raw_addr = upload["dataset_address"]

# 2. Denoise with disynthon collapse
denoise_result = denoiser.run(
    dataset_address=raw_addr,
    output_key="denoised_disynthon",
    strategy="non_unified",
    add_hit_labels=True,
    hit_percentile=90.0,
    use_disynthon_pairs=True,
    smiles_cols=["smiles_a", "smiles_b", "smiles_c"],
    aggregate_operation="sum",    # sum counts across trisynthons sharing a pair
    min_count_threshold=0,        # raise to filter low-count disynthons
)
denoised_addr = denoise_result["denoised_dataset_address"]

# 3–6. Same two-model pipeline as Workflow B,
#       but use dataset_column="disynthons" in Featurize
feat_target = featurizer.run(
    dataset_address=denoised_addr,
    featurizer="ecfp",
    output="del_di_ecfp_target",
    dataset_column="disynthons",   # ← disynthons column, not smiles
    label_column="target_hits",
)
# ... continue identically to Workflow B from here
```

### ⚠️ Warning: Disynthon + Unified

Do not combine `use_disynthon_pairs=True` with `strategy="unified"`. Count aggregation during disynthon collapse increases discrepancy between replicates, which distorts the Poisson CI calculation and produces worse enrichment scores. Always pair disynthon collapse with `strategy="non_unified"`.

---

## `DelDenoise` Parameter Reference

```python
denoiser.run(
    dataset_address,          # str     required  Datastore address of the input CSV
    output_key,               # str     required  Name for the output CSV in the datastore
    strategy,                 # str     "unified" "unified" (Poisson CI) or "non_unified" (z-score)
    control_cols,             # list    see below Control replicate column names
    target_cols,              # list    see below Target replicate column names
    add_hit_labels,           # bool    False     Append binary hit label column(s)
    hit_percentile,           # float   90.0      Percentile cutoff; top (100-N)% are hits
    alpha,                    # float   0.05      Poisson CI significance level (unified only)
    drop_duplicates,          # bool    True      Drop duplicate SMILES rows before scoring
    use_disynthon_pairs,      # bool    False     Collapse trisynthon rows into disynthon pairs
    smiles_cols,              # list    see below Three synthon SMILES columns (disynthon mode)
    aggregate_operation,      # str     "sum"     "sum" or "mean" for disynthon count aggregation
    min_count_threshold,      # int     0         Drop disynthon rows with total count below this
)
```

| Parameter | Type | Default | Notes |
|---|---|---|---|
| `dataset_address` | str | required | Must be a `deepchem://` address |
| `output_key` | str | required | `.csv` extension auto-appended if missing |
| `strategy` | str | `"unified"` | `"unified"` or `"non_unified"` |
| `control_cols` | list[str] | `["seq_matrix_1","seq_matrix_2","seq_matrix_3"]` | Background replicate columns |
| `target_cols` | list[str] | `["seq_target_1","seq_target_2","seq_target_3"]` | Target replicate columns |
| `add_hit_labels` | bool | `False` | Adds `hits` (unified) or `target_hits`+`control_hits` (non-unified) |
| `hit_percentile` | float | `90.0` | 0–100; rows above this percentile are labeled 1 |
| `alpha` | float | `0.05` | 95% CI when 0.05; unified strategy only |
| `drop_duplicates` | bool | `True` | Deduplicates on the SMILES column before scoring |
| `use_disynthon_pairs` | bool | `False` | Enables trisynthon→disynthon collapse; requires `smiles_cols` |
| `smiles_cols` | list[str] | `["smiles_a","smiles_b","smiles_c"]` | Must be exactly 3 columns |
| `aggregate_operation` | str | `"sum"` | How to combine counts for duplicate disynthon pairs |
| `min_count_threshold` | int | `0` | Filters low-evidence disynthon pairs after collapse |

**Return key:** `denoised_dataset_address` (str)

---

## Output Columns Reference

Columns added to the output CSV by `del_denoise`, depending on `strategy` and `add_hit_labels`:

| Strategy | `add_hit_labels` | Columns added |
|---|---|---|
| `unified` | `False` | `Poisson_Enrichment` — ratio of target lower CI to control upper CI; > 1 means enriched |
| `unified` | `True` | `Poisson_Enrichment`, `hits` (0/1) |
| `non_unified` | `False` | `seq_target_sum`, `seq_control_sum`, `Target_Enrichment_Score`, `Control_Enrichment_Score` |
| `non_unified` | `True` | all above + `target_hits` (0/1), `control_hits` (0/1) |

All input columns are preserved. When `use_disynthon_pairs=True`, the individual synthon SMILES columns (`smiles_a`, `smiles_b`, `smiles_c`) are replaced by a single `disynthons` column containing the RDKit-combined canonical SMILES of each pair. Use `dataset_column="disynthons"` in the downstream `Featurize` call.
