# DeepChem Server Environments

This directory contains conda environment definitions for running DeepChem Server.

## Environment Files

### core_environment.yml

The main environment for running DeepChem Server with CPU support.

**Key dependencies:**
- Python 3.11
- DeepChem (pre-release)
- PyTorch 2.1.0
- RDKit 2022.9.5
- OpenFE 1.1.0 (Free Energy calculations)
- OpenMM 8.0 (Molecular dynamics)
- VINA (Molecular docking)
- PDBFixer (Protein structure preparation)

**Installation:**

```bash
# Using mamba (recommended)
mamba env create -f core_environment.yml

# Using conda
conda env create -f core_environment.yml

# Activate
conda activate deepchem-server-env
```

### gpu_environment.yml

Environment with GPU support via CUDA toolkit.

**Additional dependencies:**
- NVIDIA CUDA Toolkit

**Installation:**

```bash
mamba env create -f gpu_environment.yml
conda activate deepchem-server-env
```

## Usage

After activating an environment, install the DeepChem Server package:

```bash
# From the repository root
pip install -e .

# Or install specific sub-packages
cd deepchem_server/core/primitives
pip install -e .
```

## Updating Environments

To update an existing environment:

```bash
mamba env update -f core_environment.yml --prune
```

To export the current environment:

```bash
conda env export > environment.yml
```

## Notes

- The `--pre` flag in pip dependencies enables installation of pre-release versions of DeepChem
- GPU environment requires NVIDIA drivers and compatible CUDA toolkit
- Some dependencies (vina, openfe, openmm, pdbfixer) are only available via conda-forge

