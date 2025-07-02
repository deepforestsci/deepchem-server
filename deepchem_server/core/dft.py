import os
import re
import ast
import tempfile
from typing import Optional, Dict, List, Iterable, Union, Set
from deepchem_server.core import config
from deepchem_server.core.cards import DataCard
import deepchem as dc
import pandas as pd
import multiprocessing as mp
from rdkit import Chem
from deepchem_server.core.address import DeepchemAddress
from gpaw import GPAW
from gpaw.eigensolvers import RMMDIIS
import ase
from ase.io import lammpsdata
import time


def calc_potential_energy(lammps_data_path: str,
                          ncores: int,
                          atomic_style: str = "atomic",
                          mode: str = "fd",
                          xc: str = "PBE",
                          convergence: dict = {}) -> ase.atoms.Atoms:
    """
    Calculate potential energy using GPAW.
    
    Parameters
    ----------
    lammps_data_path : str
        Path to the LAMMPS data file.
    ncores : int
        Number of CPU cores to use for the calculation.
    mode : str
        Calculation mode, either 'fd', 'pw', or 'lcao'.
    xc : str, optional
        Exchange-correlation functional to use, by default "pbe".
    convergence : dict, optional
        Convergence parameters for GPAW, by default {}.
    
    Returns
    -------
    None
    """
    if convergence is None or convergence == {}:
        convergence = {
            'energy': 0.0005,  # eV / electron
            'density': 1.0e-4,  # electrons / electron
            'eigenstates': 4.0e-6,  # eV^2 / electron
            'bands': 'occupied'
        }

    molecule = lammpsdata.read_lammps_data(lammps_data_path,
                                           atomic_style=atomic_style)
    # Enable periodic boundary conditions
    molecule.set_pbc(True)

    calc = GPAW(mode=mode,
                xc=xc,
                parallel={'domain': ncores},
                eigensolver=RMMDIIS(),
                h=0.18,
                convergence=convergence)

    molecule.calc = calc
    start_time = time.time()
    potential_energy = molecule.get_potential_energy()
    end_time = time.time()

    print(
        f"Potential energy calculation took {end_time - start_time:.2f} seconds"
    )
    print(f"Potential energy: {potential_energy:.6f} eV")

    return molecule


def calculate_potential_energy(lammps_file_path: str,
                               ouput_key: str,
                               atomic_style: str = "atomic",
                               mode: str = "fd",
                               xc: str = "PBE",
                               ncores: Optional[int] = None,
                               convergence: dict = {}):
    """
    Calculates the potential energy of a LAMMPS data file using GPAW and stores the updated calculator object.
    
    Parameters
    ----------
    lammps_file_path : str
        Path to the LAMMPS data file.
    ouput_key : str
        Key to store the potential energy in the calculator's data.
    atomic_style : str, optional    
        Atomic style for reading the LAMMPS data file, by default "atomic".
    mode : str, optional
        Calculation mode, either 'fd', 'pw', or 'lcao', by default "fd".
    xc : str, optional
        Exchange-correlation functional to use, by default "PBE".
    ncores : int, optional
        Number of CPU cores to use for the calculation. If None, uses only 1 core
    convergence : dict, optional
        Convergence parameters for GPAW, by default {}.
    """
    if not lammps_file_path.endswith('.data'):
        raise ValueError("LAMMPS data file must have a .data extension")

    datastore = config.get_datastore()
    tempdir = tempfile.TemporaryDirectory()
    basedir = os.path.join(tempdir.name)

    if ncores is None:
        nproc = os.cpu_count()
    else:
        nproc = ncores

    # TODO
    
    return None