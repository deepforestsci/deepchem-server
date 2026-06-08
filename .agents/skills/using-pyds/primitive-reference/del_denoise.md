# DEL Denoise

**pyds class:** `DelDenoise`
**Server primitive:** `del_denoise()` in `deepchem_server/core/primitives/del_denoising.py`
**Endpoint:** `POST /primitive/del/denoise`

## Description

Scores DNA-Encoded Library (DEL) screening data to identify enriched compounds above background noise. Accepts a CSV with per-compound sequencing counts from control (matrix/background) and target replicates, and outputs an enrichment-scored CSV.

Two scoring strategies are available:
- **Unified (Poisson CI)** — computes Poisson confidence intervals across all replicates simultaneously and derives an enrichment ratio. Statistically rigorous; recommended default.
- **Non-unified (z-score)** — sums replicate counts and computes a per-row z-score. Simpler; faster; less conservative.

Optionally collapses three-fragment trisynthon rows into all pairwise disynthon combinations before scoring.

## Inputs

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `dataset_address` | `str` | Yes | — | Datastore address of the input DEL screening CSV. |
| `output_key` | `str` | Yes | — | Short artifact name for the output scored CSV. `.csv` extension auto-appended if missing. |
| `strategy` | `str` | No | `"unified"` | Scoring strategy. `"unified"` (Poisson CI) or `"non_unified"` (z-score). |
| `control_cols` | `list[str]` | No | `["seq_matrix_1", "seq_matrix_2", "seq_matrix_3"]` | Column names holding background/control sequencing counts. |
| `target_cols` | `list[str]` | No | `["seq_target_1", "seq_target_2", "seq_target_3"]` | Column names holding target sequencing counts. |
| `add_hit_labels` | `bool` | No | `False` | If `True`, adds a binary hit column marking rows above `hit_percentile` as `1`. |
| `hit_percentile` | `float` | No | `90.0` | Percentile cutoff (0–100) for hit designation. Rows with enrichment score above this percentile are marked as hits. |
| `alpha` | `float` | No | `0.05` | Confidence level for Poisson CI calculation (unified strategy only). `0.05` = 95% CI. |
| `drop_duplicates` | `bool` | No | `True` | Remove duplicate SMILES rows before scoring. |
| `use_disynthon_pairs` | `bool` | No | `False` | If `True`, collapse trisynthon rows (3-fragment DEL) into all pairwise disynthon combinations (AB, AC, BC). |
| `smiles_cols` | `list[str]` | No | `["smiles_a", "smiles_b", "smiles_c"]` | Three SMILES column names for fragments A, B, C. Required when `use_disynthon_pairs=True`. |
| `aggregate_operation` | `str` | No | `"sum"` | How to aggregate counts when collapsing to disynthons. `"sum"` or `"mean"`. |
| `min_count_threshold` | `int` | No | `0` | Drop disynthon pairs with total count below this value (after aggregation). |

## Output

**Return key:** `"denoised_dataset_address"`
**Return type:** `str`

The returned address points to an **enrichment-scored CSV file** in the datastore.

### Output CSV columns

**Unified strategy:**

| Column | Description |
|--------|-------------|
| All input columns | Preserved from input CSV. |
| `Poisson_Enrichment` | Enrichment ratio: target lower CI bound / control upper CI bound. Values > 1.0 indicate enrichment above background. |
| `hits` | Binary (0/1) hit label. Only present if `add_hit_labels=True`. |

**Non-unified strategy:**

| Column | Description |
|--------|-------------|
| All input columns | Preserved from input CSV. |
| `seq_target_sum` | Sum of all target replicate counts. |
| `seq_control_sum` | Sum of all control replicate counts. |
| `Target_Enrichment_Score` | Z-score for target enrichment vs. uniform expectation. |
| `Control_Enrichment_Score` | Z-score for control enrichment vs. uniform expectation. |
| `target_hits`, `control_hits` | Binary hit labels. Only present if `add_hit_labels=True`. |

**Disynthon mode** (when `use_disynthon_pairs=True`):
- Original fragment SMILES columns replaced by `disynthons` column (RDKit-merged canonical SMILES of the pair).
- Count columns aggregated per disynthon pair.
- SMILES pairs that RDKit cannot parse are silently dropped.

### Metadata card fields
The data card attached to the output includes: `strategy`, `control_cols`, `target_cols`, `n_input_rows`, `n_output_rows`, `add_hit_labels`, `hit_percentile`, `alpha`, `use_disynthon_pairs`, `n_disynthons`, `aggregate_operation`, `min_count_threshold`.
