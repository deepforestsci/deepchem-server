"""
Setup script for deepchem-core-primitives package.

This package uses conda/mamba for environment management due to dependencies
on conda-forge packages (pdbfixer, vina, openfe, openmm) that are not available
via pip.

Installation:
    1. Create conda environment: mamba env create -f environment.yml
    2. Activate environment: conda activate deepchem-core-primitives
    3. Install package: pip install -e .

For development:
    pip install -e ".[dev]"
"""

from setuptools import find_packages, setup


setup(
    name="deepchem-core-primitives",
    version="0.1.0",
    description="DeepChem Server primitives - featurization, training, evaluation, inference, docking, FEP",
    author="DeepChem Team",
    python_requires=">=3.11",
    packages=find_packages(exclude=["tests", "tests.*"]),
    package_dir={"": "."},
    install_requires=[
        # Note: Most dependencies are managed via conda environment.yml
        # These are only the pip-installable dependencies
        "pandas>=2.0.0",
        "numpy>=1.24.0,<2.0",
        "scikit-learn>=1.3.2,<1.6",
        "pillow>=10.0.0",
        "deepchem",
        "torch>=2.0.0",
        "transformers>=4.30.0",
        "importlib_resources",
        # Internal dependency
        "deepchem-core-common",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
        ],
        "fep": [
            # FEP-specific dependencies (managed by conda)
            # openfe, openmm, pdbfixer are conda-only
        ],
        "docking": [
            # Docking-specific dependencies (managed by conda)
            # vina is conda-only
        ],
    },
    entry_points={
        "console_scripts": [
            # Add CLI commands if needed
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
