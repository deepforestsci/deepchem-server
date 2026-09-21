# Docking

**pyds class:** `Docking`
**Server primitive:** `generate_pose()` in `deepchem_server/core/primitives/docking.py`
**Endpoint:** `POST /primitive/generate_pose`

## Description

Performs automated molecular docking using AutoDock VINA via DeepChem's `VinaPoseGenerator`. Given a protein structure and a ligand, generates multiple predicted binding poses ranked by binding affinity (kcal/mol). Outputs a results summary JSON containing addresses to per-mode PDB complexes and a binding scores JSON.

**Dependencies:** Requires AutoDock VINA and RDKit to be installed in the server environment.

## Inputs

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `protein_address` | `str` | Yes | — | Datastore address of the protein PDB file. Must be a prepared protein structure (hydrogens present, no clashing atoms). Use `PdbClean` output for best results. |
| `ligand_address` | `str` | Yes | — | Datastore address of the ligand file. Accepts `.pdb` or `.sdf` format — format detected from file extension. |
| `output` | `str` | Yes | — | Prefix for all generated output files (e.g., `"dock_result"` produces `dock_result_mode_1.pdb`, `dock_result_scores.json`, `dock_result_results.json`). |
| `exhaustiveness` | `int` | No | `10` | VINA search exhaustiveness. Higher values increase search thoroughness and computation time. Typical range: 1–32. |
| `num_modes` | `int` | No | `9` | Number of binding poses to generate. VINA may produce fewer if convergence issues occur; remaining slots are filled with the last valid pose's score. |
| `save_pdbqt` | `bool` | No | `False` | Whether to also save PDBQT-format ligand pose files (preserve atom charges and partial bonds for downstream VINA analysis). PDB complexes are always generated. |

## Output

**Return key:** `"docking_results_address"`
**Return type:** `str`

The returned address points to a **results summary JSON file** (`{output}_results.json`).

### Results summary JSON structure

```json
{
  "docking_method": "VINA",
  "exhaustiveness": 10,
  "complex_addresses": {
    "mode 1": "deepchem://profile/project/dock_result_mode_1.pdb",
    "mode 2": "deepchem://profile/project/dock_result_mode_2.pdb"
  },
  "scores_address": "deepchem://profile/project/dock_result_scores.json",
  "pdbqt_addresses": {
    "mode 1": "deepchem://profile/project/dock_result_mode_1.pdbqt"
  },
  "message": "VINA docking completed successfully"
}
```

`pdbqt_addresses` is only present if `save_pdbqt=True`.

### Scores JSON structure (at `scores_address`)

```json
{
  "mode 1": {"affinity (kcal/mol)": -8.4},
  "mode 2": {"affinity (kcal/mol)": -7.9}
}
```

Binding affinity in kcal/mol — **more negative = stronger predicted binding**.

### Files created in datastore

| File | Description |
|------|-------------|
| `{output}_mode_{i}.pdb` | Protein-ligand complex PDB for pose i (RDKit `CombineMols` output). For visualization. |
| `{output}_scores.json` | Binding affinity for each mode in kcal/mol. |
| `{output}_mode_{i}.pdbqt` | PDBQT ligand pose for mode i (only if `save_pdbqt=True`). For VINA re-scoring or downstream analysis. |
| `{output}_results.json` | Master summary file. This is what the return address points to. |

### Notes
- PDB complexes are created via RDKit's `CombineMols` — they are suitable for visualization but lack the full VINA charge/bond model. For detailed pose analysis, use PDBQT files.
- If fewer poses are generated than `num_modes`, the last valid affinity score is repeated for the remaining mode slots in the scores JSON.
- Grid box for docking is auto-generated around the ligand's initial position. Ensure the ligand is placed in or near the binding site in the input file.
