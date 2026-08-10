# Partition

**pyds class:** `Partition`
**Server primitive:** `partition()` in `deepchem_server/core/primitives/partition.py`
**Endpoint:** `POST /primitive/partition`

## Description

Divides a dataset into N roughly equal, non-overlapping partitions. Primarily used to split large datasets before parallel featurization or distributed training. Each partition is uploaded to the datastore as an independent dataset with metadata linking it to its parent.

## Inputs

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `dataset_address` | `str` | Yes | — | Datastore address of the dataset to partition. Supports `DiskDataset` and CSV/DataFrame. |
| `n_partition` | `int` | No | `4` | Number of partitions to create. Must be > 0. |
| `shuffle` | `bool` | No | `False` | Whether to randomly shuffle samples before partitioning. **Only supported for `DiskDataset`** — raises `ValueError` if set to `True` on a CSV dataset. |

## Output

**Return key:** `"partitioned_dataset_addresses"`
**Return type:** `list[str]` — length equals `n_partition`

```python
[
  "deepchem://profile/project/<key>_partition0",
  "deepchem://profile/project/<key>_partition1",
  ...
]
```

Each partition's data card includes:
- `n_partition` — total number of partitions
- `partition_id` — zero-based index of this partition
- `parent` — datastore address of the original dataset

The parent dataset's card is also updated with `n_partition` to mark it as partitioned.

### Partition sizing
- **DiskDataset:** Partitions are approximately equal. Due to integer division, later partitions may have slightly fewer samples.
- **CSV:** Partitions are streamed in chunks of `min(100_000, rows_per_partition)` rows for memory efficiency. Boundaries may not align perfectly with target partition size.

### Featurizer metadata
- Preserved for `DiskDataset` partitions.
- Not carried over for CSV partitions.

### Common use case
Partition a large raw CSV → featurize each partition independently (in parallel or sequentially) → merge the resulting DiskDatasets.
