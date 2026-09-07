# PDB Clean

**pyds class:** `PdbClean` (import from `pyds.primitives`, not top-level `pyds`)
**Server primitive:** `pdb_clean()` in `deepchem_server/core/primitives/proteome_scan/`
**Endpoint:** `POST /primitive/pdb_clean`

## Description

Prepares protein structures for computational workflows (docking, FEP) by:
1. Querying the PDBe API for all experimental structures matching a UniProt entry.
2. Selecting a non-overlapping, high-resolution set of structures that maximizes sequence coverage.
3. Cleaning each selected PDB: removing non-target chains, stripping all heterogens and water, adding hydrogens at pH 7.0.
4. Uploading the cleaned PDBs and a summary JSON to the datastore.

**Dependencies:** Requires PDBFixer, OpenMM, and BioPython in the server environment. Requires network access to the PDBe EMBL-EBI API.

## Inputs

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `gene_name` | `str` | Yes | — | Human gene name (HUGO symbol), e.g., `"MAP2K1"`, `"MB"`. Used for cache file naming. |
| `entry_id` | `str` | Yes | — | UniProt accession ID, e.g., `"P02144"`. Used to query PDBe for all associated experimental structures. |
| `scan_id` | `str` | Yes | — | Unique identifier for this scan run, e.g., `"scan_001"`. All cache artifacts are namespaced under `PROTEOMESCAN_CACHE_ROOT/<scan_id>/<gene_name>/`. |
| `output` | `str` | Yes | — | Directory prefix for all output files in the datastore (e.g., `"mb_results"`). All files are stored under `deepchem://profile/project/<output>/`. |
| `min_res_val` | `float` | No | `2.5` | Resolution threshold in Ångströms for the high-resolution candidate set. Only structures with resolution ≤ this value are considered in the first selection pass. |

## Output

**Return key:** `"pdb_clean_results_address"`
**Return type:** `str`

The returned address points to a **JSON summary file** in the datastore.

### Summary JSON structure

```json
{
  "gene_name": "MB",
  "entry_id": "P02144",
  "scan_id": "scan_001",
  "min_res_val": 2.5,
  "selected_pdb_ids": ["1MBN", "2MBN"],
  "cleaned_pdb_addresses": {
    "1MBN": "deepchem://profile/project/mb_results/cleaned_g_MB_p_1MBN.pdb",
    "2MBN": "deepchem://profile/project/mb_results/cleaned_g_MB_p_2MBN.pdb"
  },
  "pdbs_metadata_csv_address": "deepchem://profile/project/mb_results/MB_pdbs.csv"
}
```

### Files created in datastore

| File pattern | Description |
|---|---|
| `{output}/cleaned_g_{gene_name}_p_{pdb_id}.pdb` | Cleaned PDB file for each selected structure. Heterogens and water removed; hydrogens added at pH 7.0. |
| `{output}/{gene_name}_pdbs.csv` | Metadata table of selected PDBs (PDB ID, resolution, chain coverage). |
| `{output}/{gene_name}_pdb_clean_summary.json` | Summary JSON — the address this primitive returns. |

### Metadata CSV columns

| Column | Description |
|--------|-------------|
| `id` | PDB ID (uppercase, e.g., `"1MBN"`) |
| `resolution` | X-ray resolution in Ångströms |
| `chains` | Chain IDs and residue range selected (e.g., `"A/B=10-290"`) |
| `chain_start`, `chain_end` | Residue index boundaries |
| `coverage` | Sequence span covered (`chain_end - chain_start`) |

### PDB selection algorithm
- Structures are selected greedily to maximize sequence coverage with minimal overlap.
- Two candidate sets are built: structures ≤ `min_res_val` Å, and structures ≤ the resolution of the highest-coverage structure.
- Within overlapping ranges, a structure replaces a prior selection if it has similar resolution (< 0.5 Å difference) and greater coverage.
- The union of both candidate sets forms the final selection.

### Cleaning steps applied to each PDB
1. Remove all chains not associated with the target protein.
2. Remove all heterogens (ligands, ions) and water (`keepWater=False`).
3. Add missing hydrogens at physiological pH 7.0 using PDBFixer/OpenMM.

### Notes
- Raw (uncleaned) PDB files are cached locally at `PROTEOMESCAN_CACHE_ROOT/<scan_id>/<gene_name>/`. Set this environment variable before running.
- If a PDB fails to clean (PDBFixer error), it is dropped and the next best structure is used. If all PDBs fail, a `RuntimeError` is raised.
- **No ligands are retained** in the output — all heterogens are stripped. If a co-crystallized ligand is needed as a docking reference, it must be extracted before calling this primitive.
- The cleaned PDB addresses (from `cleaned_pdb_addresses`) are the inputs to `Docking` or `run_rbfe`.
