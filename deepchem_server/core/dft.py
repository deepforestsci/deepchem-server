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
                          ncores: Optional[int] = None,
                          atomic_style: str = "atomic",
                          mode: str = "fd",
                          xc: str = "PBE",
                          convergence: Optional[dict] = None) -> ase.atoms.Atoms:
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

    return calc


def calculate_potential_energy(datafile_address: str,
                               output_key: str,
                               atomic_style: str = "atomic",
                               mode: str = "fd",
                               xc: str = "PBE",
                               ncores: Optional[int] = None,
                               convergence: Optional[dict] = None):
    """
    Calculates the potential energy of a LAMMPS data file using GPAW and stores the updated calculator object.
    
    Parameters
    ----------
    lammps_file_path : str
        Path to the LAMMPS data file.
    output_key : str
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
    if isinstance(datafile_address, str):
        datafile_address = ast.literal_eval(datafile_address)
    
    if isinstance(output_key, str):
        output_key = ast.literal_eval(output_key)
    
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore is not configured. Please set up the datastore.")

    datafile_path = datastore.get_path(datafile_address)
    calc_obj = calc_potential_energy(
        lammps_data_path=datafile_path,
        ncores=ncores,
        atomic_style=atomic_style,
        mode=mode,
        xc=xc,
        convergence=convergence
    )
    
    description = f"ASE Atoms object with attached GPAW calculator for potential energy calculation from LAMMPS data file: {datafile_address}"
    card = DataCard(
        address='',
        file_type='gpw',
        data_type='ase.atoms.Atoms',
        description=description
    )
    tempdir = tempfile.TemporaryDirectory()
    temp_output_path = os.path.join(tempdir.name, 'temp.gpw')
    if not output_key.endswith('.gpw'):
        output_key += '.gpw'

    calc_obj.write(temp_output_path, mode='all')
    output_address = datastore.upload_data(
        DeepchemAddress.get_key(output_key), temp_output_path, card
    )
    return output_address