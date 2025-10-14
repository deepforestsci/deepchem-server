import os
import json
import tempfile
from deepchem_server.core import config
from deepchem_server.core.cards import DataCard
from deepchem_server.core.progress_logger import log_progress
from deepchem.dock.pose_generation import VinaPoseGenerator

# Import RDKit at module level
try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


def split_pdbqt_docked_ligands(docked_ligand_name: str, docking_method: str):
    """
    Split PDBQT file with multiple modes into separate files for each mode.
    """
    with open(docked_ligand_name, 'r+') as fp:
        # Read and store all lines into list
        lines = fp.readlines()

        # Remove the MODEL 'n' and ENDMDL lines when num_modes is greater than 1
        ind_MDL = [i for i, x in enumerate(lines) if "MODEL" in x]
        ind_ENDMDL = [i for i, x in enumerate(lines) if "ENDMDL" in x]
        remove_ind = ind_MDL + ind_ENDMDL
        remove_ind.sort()
        lines = [i for j, i in enumerate(lines) if j not in remove_ind]

        if docking_method == 'VINA' or docking_method == 'qVINA-W':
            # Find the indexes to split the ligands
            indx = [i for i, x in enumerate(lines) if "REMARK VINA RESULT:" in x]
            indx.append(len(lines) - 1)

        elif docking_method == 'GNINA':
            indx = [i for i, x in enumerate(lines) if "REMARK  Name =" in x]
            indx.append(len(lines) - 1)

        # Write the pdbqt docked ligands
        for i in range(len(indx)):
            if i == len(indx) - 1:
                with open("%s_ligand_docked.pdbqt" % str(i), 'w') as fp_i:
                    # Move file pointer to the beginning of a file
                    fp_i.seek(0)
                    # Truncate the file
                    fp_i.truncate()
                    # Remove lines (first and last) that breaks Pose Scoring code
                    fp_i.writelines(lines[indx[i - 1]:indx[i] + 1])
            else:
                with open(("%s_ligand_docked.pdbqt" % str(i + 1)), 'w') as fp:
                    # Move file pointer to the beginning of a file
                    fp.seek(0)
                    # Truncate the file
                    fp.truncate()
                    # Remove lines (first and last) that breaks Pose Scoring code
                    fp.writelines(lines[indx[i]:indx[i + 1]])


def generate_pose(
    protein_address: str,
    ligand_address: str,
    output: str,
    exhaustiveness: int = 10,
    num_modes: int = 9,
    save_pdbqt: bool = False,
) -> str:
    """
    Generate VINA molecular docking poses.

    Parameters
    ----------
    protein_address: str
        DeepChem address of the protein PDB file
    ligand_address: str
        DeepChem address of the ligand file (PDB or SDF)
    output: str
        Output name for the docking results
    exhaustiveness: int
        Vina exhaustiveness parameter (default: 10)
    num_modes: int
        Number of binding modes to generate (default: 9)
    save_pdbqt: bool
        Whether to save PDBQT files in addition to PDB complexes (default: False)

    Returns
    -------
    str
        DeepChem address of the docking results
    """

    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")

    if not protein_address or not ligand_address:
        raise ValueError('Protein and/or ligand input is required.')

    try:
        # Check dependencies
        if not RDKIT_AVAILABLE:
            raise ImportError("RDKit is required for docking but not installed")

        # Check if VINA is available (VinaPoseGenerator will fail if not)
        try:
            pg = VinaPoseGenerator()
        except Exception as e:
            raise ImportError(f"VINA/AutoDock VINA is required for docking but not available: {e}")

        with tempfile.TemporaryDirectory() as tmp:
            log_progress('docking', 10, f'downloading protein from {protein_address}')
            protein_path = os.path.join(tmp, 'protein.pdb')
            datastore.download_object(protein_address, protein_path)

            log_progress('docking', 20, f'downloading ligand from {ligand_address}')
            # Detect format from address and let DeepChem handle conversion
            ligand_ext = '.sdf' if ligand_address.endswith('.sdf') else '.pdb'
            ligand_path = os.path.join(tmp, f'ligand{ligand_ext}')
            datastore.download_object(ligand_address, ligand_path)

            log_progress('docking', 30, 'preparing molecules for VINA')

            log_progress('docking', 40, 'initializing VINA pose generator')
            log_progress('docking', 50, f'generating {num_modes} poses with VINA')
            # Generate poses using file paths - DeepChem handles preparation internally
            complexes, scores = pg.generate_poses(molecular_complex=(protein_path, ligand_path),
                                                  exhaustiveness=exhaustiveness,
                                                  num_modes=num_modes,
                                                  out_dir=tmp,
                                                  generate_scores=True)

            # Validate that we got valid results
            if not complexes or not scores:
                raise ValueError("No docking poses or scores generated")

            # Ensure we don't exceed available results
            actual_modes = min(num_modes, len(complexes), len(scores))
            if actual_modes == 0:
                raise ValueError("No valid docking results generated")

            log_progress('docking', 60, f'generated {actual_modes} valid poses')

            # Handle PDBQT files if requested
            pdbqt_addresses = {}
            if save_pdbqt:
                log_progress('docking', 65, 'processing PDBQT files')
                docked_ligand_name = os.path.join(tmp, "temp_ligand_docked.pdbqt")

                if os.path.exists(docked_ligand_name):
                    # Clean the PDBQT file (remove first and last lines that break pose scoring)
                    with open(docked_ligand_name, 'r+') as fp:
                        lines = fp.readlines()
                        fp.seek(0)
                        fp.truncate()
                        fp.writelines(lines[1:-1])

                    # Split PDBQT file if multiple modes
                    if actual_modes > 1:
                        # Change to temp directory for split files
                        original_cwd = os.getcwd()
                        os.chdir(tmp)
                        split_pdbqt_docked_ligands(docked_ligand_name, "VINA")
                        os.chdir(original_cwd)

                    # Upload PDBQT files for each mode
                    for i in range(actual_modes):
                        try:
                            pdbqt_filename = f"{output}_mode_{i + 1}.pdbqt"
                            pdbqt_card = DataCard(address='', file_type='pdbqt', data_type='text/plain')

                            if actual_modes == 1:
                                # Single mode: use the original file
                                pdbqt_content = open(docked_ligand_name, 'r').read()
                            else:
                                # Multiple modes: use the split files
                                split_filename = os.path.join(tmp, f"{i + 1}_ligand_docked.pdbqt")
                                if os.path.exists(split_filename):
                                    pdbqt_content = open(split_filename, 'r').read()
                                    # Clean up the split file
                                    os.remove(split_filename)
                                else:
                                    log_progress('docking', 68, f'Warning: Split PDBQT file not found for mode {i + 1}')
                                    continue

                            pdbqt_address = datastore.upload_data_from_memory(pdbqt_content, pdbqt_filename, pdbqt_card)

                            if pdbqt_address:
                                pdbqt_addresses['mode %s' % (i + 1)] = pdbqt_address
                                log_progress('docking', 67, f'saved PDBQT for mode {i + 1}')
                        except Exception as e:
                            log_progress('docking', 68, f'failed to save PDBQT for mode {i + 1}: {e}')
                else:
                    log_progress('docking', 66, 'Warning: PDBQT file not found in temp directory')

            log_progress('docking', 70, 'preparing results')
            # Format scores: always include requested mode keys; pad with last available score if needed
            scores_formatted = {}
            complex_addresses = {}
            modes_to_report = max(actual_modes, num_modes)

            for i in range(modes_to_report):
                idx = min(i, actual_modes - 1)
                scores_formatted['mode %s' % (i + 1)] = {'affinity (kcal/mol)': float(scores[idx])}

                # Save complex PDB file for each pose
                if idx < len(complexes) and complexes[idx] is not None:
                    try:
                        # Combine protein and ligand from the complex
                        complex_mol = Chem.CombineMols(complexes[idx][0], complexes[idx][1])

                        # Create complex file content
                        complex_content = Chem.MolToPDBBlock(complex_mol)

                        # Upload complex file
                        complex_filename = f"{output}_mode_{i + 1}.pdb"
                        complex_card = DataCard(address='', file_type='pdb', data_type='text/plain')
                        complex_address = datastore.upload_data_from_memory(complex_content, complex_filename,
                                                                            complex_card)

                        if complex_address:
                            complex_addresses['mode %s' % (i + 1)] = complex_address
                            log_progress('docking', 75, f'saved complex for mode {i + 1}')
                    except Exception as e:
                        log_progress('docking', 76, f'failed to save complex for mode {i + 1}: {e}')

            # Upload a standalone scores JSON and capture its datastore address
            try:
                scores_card = DataCard(address='', file_type='json', data_type='json')
                scores_json_str = json.dumps(scores_formatted)
                scores_address = datastore.upload_data_from_memory(scores_json_str, f"{output}_scores.json",
                                                                   scores_card)
            except Exception as e:
                scores_address = None
                log_progress('docking', 72, f'failed to upload scores JSON: {e}')

            results = {
                'docking_method': 'VINA',
                'exhaustiveness': exhaustiveness,
                'complex_addresses': complex_addresses,
                'scores_address': scores_address,
                'message': 'VINA docking completed successfully',
            }

            # Add PDBQT addresses if they were generated
            if pdbqt_addresses:
                results['pdbqt_addresses'] = pdbqt_addresses

            log_progress('docking', 90, 'uploading results summary')
            # Upload results summary: file is JSON, logical data type is 'json'
            card = DataCard(address='', file_type='json', data_type='json')
            results_json = json.dumps(results)
            result_address = datastore.upload_data_from_memory(results_json, f"{output}_results.json", card)

            if result_address is None:
                raise ValueError("Failed to upload docking results to datastore")

            log_progress('docking', 100, 'VINA docking completed successfully')
            return result_address

    except Exception as e:
        raise Exception(f'VINA docking failed: {str(e)}')
