# RBFE (Relative Binding Free Energy)

**pyds class:** Not yet exposed via pyds client (server-side only)
**Server primitive:** `run_rbfe()` and `collate_rbfe_results()` in `deepchem_server/core/primitives/fep/rbfe/`
**Endpoints:** `POST /primitive/run_rbfe`, `POST /primitive/collate_rbfe_results`

## Description

Computes **relative binding free energy (RBFE)** differences between pairs of ligands using alchemical free energy perturbation (FEP) with the `RelativeHybridTopologyProtocol` from OpenFE. The workflow has two stages:

1. **`run_rbfe`** — generates a perturbation network between ligands, sets up and runs one or more alchemical MD simulations (edges), and writes per-edge results (ΔΔG, uncertainties) to JSON files.
2. **`collate_rbfe_results`** — reads all per-edge JSON results, propagates ΔΔG values across the network from a reference ligand, and writes a ranked CSV of absolute ΔG estimates.

**Dependencies:** Requires OpenFE, OpenMM, and RDKit in the server environment.

---

## `run_rbfe` — Run RBFE Simulations

### Inputs

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `ligands_sdf_address` | `str` | Yes | — | Datastore address of the ligands SDF file. Must contain valid 3D coordinates and explicit hydrogens. Ligand names in the SDF become component names in results. |
| `cleaned_protein_pdb_address` | `str` | Yes | — | Datastore address of the cleaned protein PDB. Use `PdbClean` output. |
| `network_type` | `str` | Yes | — | Perturbation network topology. See valid values below. |
| `scorer_type` | `str` or `None` | Yes | — | Similarity scorer for mapping ligand pairs. `"LOMAP"` or `None` (only valid for `RADIAL` networks). |
| `solvent_json` | `str` | Yes | — | JSON string defining solvent composition. See format below. |
| `overridden_rbfe_settings` | `str` | Yes | — | JSON string with simulation parameter overrides. Can be `"{}"` to use defaults. |
| `dry_run` | `bool` | No | `False` | If `True`, validates setup without running MD simulations. Result JSONs will have `is_dry_run: true`. |
| `radial_network_central_ligand` | `str` or `None` | No | `None` | Name of the central (hub) ligand for `RADIAL` networks. If `None`, auto-selected as the ligand with the highest average LOMAP similarity to all others. |
| `output_key` | `str` or `None` | No | `None` | Directory prefix for result JSON files. Defaults to the parent directory of `ligands_sdf_address`. |

### `solvent_json` format

```json
{
  "positive_ion": "Na",
  "negative_ion": "Cl",
  "ion_concentration": "0.15 mol/L",
  "neutralize": true
}
```

### `overridden_rbfe_settings` common keys

```json
{
  "simulation_settings.equilibration_length": "2 ns",
  "simulation_settings.production_length": "10 ns",
  "protocol_repeats": 3
}
```

Default production length is 10 ns. Set shorter for testing.

### `network_type` valid values

| Value | Edges | Scorer required | Use when |
|-------|-------|----------------|----------|
| `"RADIAL"` | N−1 | No (LOMAP used for auto-selection only) | Hub-and-spoke: one central reference ligand, fast, fewer simulations. |
| `"MINIMAL_SPANNING"` | N−1 | Yes (`"LOMAP"`) | Minimum spanning tree over similarity graph. Balanced coverage with few edges. |
| `"MAXIMAL"` | N×(N−1)/2 | Yes (`"LOMAP"`) | All pairwise edges. Highest redundancy; enables cycle-closure error estimation. |

### Output

**Return key:** `"relative_binding_free_energy_results_address"`
**Return type:** `list[str]` — one address per edge simulated

Each address points to a **per-edge result JSON file** named `rbfe-{componentA}-{componentB}.json`.

### Per-edge JSON structure

```json
{
  "componentA_name": "benzene",
  "componentB_name": "toluene",
  "is_dry_run": false,
  "complex_dG": "-1.45 kilocalorie / mole",
  "complex_dG_uncertainty": "0.12 kilocalorie / mole",
  "solvent_dG": "-0.32 kilocalorie / mole",
  "solvent_dG_uncertainty": "0.09 kilocalorie / mole",
  "edge_ddG": "-1.13 kilocalorie / mole",
  "edge_ddG_uncertainty": "0.15 kilocalorie / mole"
}
```

- `complex_dG` — ΔG of the alchemical transformation in the protein-bound state.
- `solvent_dG` — ΔG of the transformation in solvent only (unbound).
- `edge_ddG` = `complex_dG − solvent_dG` — the **relative binding free energy difference** between ligand A and B.
- Uncertainties are propagated as root-sum-square: `√(complex_unc² + solvent_unc²)`.

### AWS Batch parallel execution
If `AWS_BATCH_JOB_ARRAY_INDEX` is set in the environment, only the edge at that array index is executed. This enables each edge to run as a separate AWS Batch array job for massively parallel RBFE campaigns.

---

## `collate_rbfe_results` — Aggregate Edge Results

### Inputs

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `result_file_deepchem_addresses` | `list[str]` | Yes | List of all per-edge result JSON addresses from `run_rbfe`. |
| `reference_ligand` | `str` or `None` | No | Name of a ligand with known experimental ΔG. Used to anchor relative values to absolute scale. If `None`, first ligand in results is used with ΔG = 0 kcal/mol. |
| `reference_ligand_dg_value` | `str` or `None` | No | Experimental ΔG of the reference ligand, e.g., `"-9.0 kcal/mol"`. Required if `reference_ligand` is set. |
| `reference_ligand_dg_value_uncertainty` | `str` or `None` | No | Uncertainty in the reference ΔG, e.g., `"0.5 kcal/mol"`. Required if `reference_ligand` is set. |
| `output_file_name` | `str` | Yes | Filename for the output CSV. `.csv` extension auto-appended if missing. |

### Output

**Return key:** `"collate_relative_binding_free_energy_results_address"`
**Return type:** `str`

The returned address points to a **ranked CSV file**.

### Output CSV structure

```
Ligand,DG Value,Uncertainty
toluene,-10.25 kilocalorie / mole,0.38 kilocalorie / mole
benzene,-9.12 kilocalorie / mole,0.35 kilocalorie / mole
xylene,-8.87 kilocalorie / mole,0.42 kilocalorie / mole
```

Rows are sorted by ΔG value ascending (most negative = tightest binder at the top).

If no reference ligand with an experimental value is provided, ΔG values are relative (anchored to 0 kcal/mol) and indicate **rank order**, not absolute binding affinity.

### ΔG propagation
- Values are propagated iteratively through the perturbation network graph from the reference ligand.
- Uncertainties accumulate additively along each path (assumes independent errors).
- Ligands not connected to the reference through valid edges will have no ΔG assigned.
