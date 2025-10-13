import os
import json
import tempfile
from deepchem_server.core import config
from deepchem_server.core.cards import DataCard
from deepchem_server.core.progress_logger import log_progress
from deepchem.dock.pose_generation import VinaPoseGenerator


def generate_pose(
    protein_address: str,
    ligand_address: str,
    output: str,
    exhaustiveness: int = 10,
    num_modes: int = 9,
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
        tempdir = tempfile.TemporaryDirectory()

        log_progress('docking', 10, f'downloading protein from {protein_address}')
        protein_path = os.path.join(tempdir.name, 'protein.pdb')
        datastore.download_object(protein_address, protein_path)

        log_progress('docking', 20, f'downloading ligand from {ligand_address}')
        # Detect format from address and let DeepChem handle conversion
        ligand_ext = '.sdf' if ligand_address.endswith('.sdf') else '.pdb'
        ligand_path = os.path.join(tempdir.name, f'ligand{ligand_ext}')
        datastore.download_object(ligand_address, ligand_path)

        log_progress('docking', 30, 'preparing molecules for VINA')

        log_progress('docking', 40, 'initializing VINA pose generator')
        pg = VinaPoseGenerator()

        with tempdir as tmp:
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

            log_progress('docking', 70, 'preparing results')
            # Format scores: always include requested mode keys; pad with last available score if needed
            scores_formatted = {}
            complex_addresses = {}
            modes_to_report = max(actual_modes, num_modes)

            for i in range(modes_to_report):
                idx = min(i, actual_modes - 1)
                scores_formatted['mode %s' % (i + 1)] = {'affinity (kcal/mol)': float(scores[idx])}
                
                # Save complex PDB file for each pose (like Chiron does)
                if idx < len(complexes) and complexes[idx] is not None:
                    try:
                        from rdkit import Chem
                        # Combine protein and ligand from the complex (same as Chiron)
                        complex_mol = Chem.CombineMols(complexes[idx][0], complexes[idx][1])
                        
                        # Create complex file content
                        complex_content = Chem.MolToPDBBlock(complex_mol)
                        
                        # Upload complex file
                        complex_filename = f"{output}_mode_{i+1}.pdb"
                        complex_card = DataCard(address='', file_type='pdb', data_type='text/plain')
                        complex_address = datastore.upload_data_from_memory(
                            complex_content, complex_filename, complex_card)
                        
                        if complex_address:
                            complex_addresses['mode %s' % (i + 1)] = complex_address
                            log_progress('docking', 75, f'saved complex for mode {i+1}')
                    except Exception as e:
                        log_progress('docking', 76, f'failed to save complex for mode {i+1}: {e}')

            results = {
                'docking_method': 'VINA',
                'num_modes': actual_modes,
                'scores': scores_formatted,
                'complexes_count': len(complexes),
                'complex_addresses': complex_addresses,
                'message': 'VINA docking completed successfully',
            }

            log_progress('docking', 90, 'uploading results summary')
            # Upload results summary: file is JSON, logical data type is 'docking results'
            card = DataCard(address='', file_type='json', data_type='docking results')
            results_json = json.dumps(results)
            result_address = datastore.upload_data_from_memory(results_json, f"{output}_results.json", card)

            if result_address is None:
                raise ValueError("Failed to upload docking results to datastore")

            log_progress('docking', 100, 'VINA docking completed successfully')
            return result_address

    except Exception as e:
        raise Exception(f'VINA docking failed: {str(e)}')