import os
import subprocess
from pathlib import Path
import MDAnalysis as mda
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from openmm.app import *
from openmm import *
from openmm.unit import *

def run_command(cmd):
    """
    Runs a shell command and returns its output, error, and exit code.
    """
    result = subprocess.run(
        cmd, 
        shell=True,           # Use True if cmd is a string
        capture_output=True,  # Capture stdout and stderr
        text=True             # Decode bytes to string
    )
    
    return {
        "stdout": result.stdout.strip(),
        "stderr": result.stderr.strip(),
        "returncode": result.returncode
    }

PATH = os.path.dirname(os.path.abspath(__file__))

# Step 1: removed water
# TODO: Yet to handle crystalline water
def remove_water(input_pdb_path, output_pdb_path):
    command = f"grep -v HOH {input_pdb_path} > {output_pdb_path}"
    output = run_command(command)
    print("remove water logs:")
    print(output)


# Step 2: isolate protein (the input protein has all the waters and other stuff removed. Also remove all the Hydrogens)
def isolate_protein(input_pdb_path, output_pdb_path):
    u = mda.Universe(input_pdb_path)
    protein = u.select_atoms("protein and not name H* 1H* 2H* 3H*").write(output_pdb_path)

# Step 3: add caps (if required)
def add_caps(input_pdb_path, output_pdb_path, is_multi_chain=False):
    capping_command = f"python {PATH}/add_caps.py -i {input_pdb_path} -o {output_pdb_path}"
    output = run_command(capping_command)
    print("add caps logs:")
    print(output)

    if is_multi_chain:
        print("adding TER between multiple chains")
        TER_command = """awk '/NME/{nme=NR} /ACE/ && nme && NR > nme {print "TER"; nme=0} {print}'""" + f" {output_pdb_path} > {output_pdb_path}"
        output = run_command(TER_command)
        print(output)

# Step 4: fix capped protien
def fix_protien(input_pdb_path, output_pdb_path):

    # clean up the original PDB file and add missing residues and heavy atoms
    fixer = PDBFixer(input_pdb_path)
    fixer.findMissingResidues()

    # only add missing residues in the middle of the chain, do not add terminal ones
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    missingResidues = dict()
    for key in keys:
        chain = chains[key[0]]
        if not (key[1] == 0 or key[1] == len(list(chain.residues()))):
            missingResidues[key] = fixer.missingResidues[key]
    fixer.missingResidues = missingResidues

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb_path, 'w'))


# Step 5: Add Hydrogen (optionally minimize)
def add_hydrogen(input_pdb_path, output_pdb_path, run_minimization_test=False, minimized_pdb_output_path=None):
    pdb = PDBFile(input_pdb_path)

    # Specify the forcefield
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)

    # Open a file in write mode and save the PDB
    with open(output_pdb_path, 'w') as output_file:
        PDBFile.writeFile(modeller.topology, modeller.positions, output_file)

    if run_minimization_test:
        from sys import stdout
        if minimized_pdb_output_path is None:
            raise ValueError("Please provide value for 'minimized_pdb_output_path'")
        modeller.addSolvent(forcefield, padding=1.0*nanometer)
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        print("Minimizing energy")
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter(minimized_pdb_output_path, 1000))
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                potentialEnergy=True, temperature=True, volume=True))
        simulation.reporters.append(StateDataReporter("md_log.txt", 100, step=True,
                potentialEnergy=True, temperature=True, volume=True))
        print("Running NVT")
        simulation.step(10000)


def prepare_protein(input_pdb_path, is_multi_chain=False, run_minimization_test=False):
    """
    """
    path = Path(input_pdb_path)
    input_pdb_name = path.stem
    root = path.parent

    output_pdb_name_1 = input_pdb_name + "_without_H2O"
    output_pdb_path_1 = os.path.join(root, output_pdb_name_1 + ".pdb")
    remove_water(input_pdb_path, output_pdb_path_1)

    output_pdb_name_2 = output_pdb_name_1 + "_protein"
    output_pdb_path_2 = os.path.join(root, output_pdb_name_2 + ".pdb")
    isolate_protein(output_pdb_path_1, output_pdb_path_2)
    
    output_pdb_name_3 = output_pdb_name_2 + "_capped"
    output_pdb_path_3 = os.path.join(root, output_pdb_name_3 + ".pdb")
    add_caps(output_pdb_path_2, output_pdb_path_3, is_multi_chain=is_multi_chain)

    output_pdb_name_4 = output_pdb_name_3 + "_fixed"
    output_pdb_path_4 = os.path.join(root, output_pdb_name_4 + ".pdb")
    fix_protien(output_pdb_path_3, output_pdb_path_4)

    output_pdb_name_5 = output_pdb_name_4 + "_with_H"
    output_pdb_path_5 = os.path.join(root, output_pdb_name_5 + ".pdb")

    if run_minimization_test:
        minimized_output_pdb_name_5 = output_pdb_name_5 + "_minimized"
        minimized_output_pdb_path_5 = os.path.join(root, minimized_output_pdb_name_5 + ".pdb")
    else:
        minimized_output_pdb_path_5 = None
    add_hydrogen(output_pdb_path_4, output_pdb_path_5, run_minimization_test=run_minimization_test, minimized_pdb_output_path=minimized_output_pdb_path_5)

    return output_pdb_path_5

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_pdb_path", type=str, required=True)
    parser.add_argument("--is_multi_chain", type=bool, default=False)
    parser.add_argument("--run_minimization_test", type=bool, default=False)

    args = parser.parse_args()
    input_pdb_path = args.input_pdb_path
    is_multi_chain = args.is_multi_chain
    run_minimization_test = args.run_minimization_test
    output_pdb = prepare_protein(input_pdb_path, is_multi_chain=is_multi_chain, run_minimization_test=run_minimization_test)
    print("prepped protein PDB: ", output_pdb)
