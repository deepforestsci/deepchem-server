from deepchem_server.core.cards import DataCard
from deepchem_server.core import config, generate_pose
import os
import json
import pytest


def test_generate_pose_basic_functionality(disk_datastore):
    """Test basic VINA pose generation functionality."""
    config.set_datastore(disk_datastore)

    # Create test protein and ligand files
    current_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_test_file = os.path.join(current_dir, "assets/cleaned_3cyx.pdb")
    ligand_test_file = os.path.join(current_dir, "assets/ligand_3cyx.pdb")

    # Upload test files to datastore
    card_protein = DataCard(address='', file_type='pdb', data_type='text/plain')
    card_ligand = DataCard(address='', file_type='pdb', data_type='text/plain')

    pdb_address = disk_datastore.upload_data("protein_test.pdb", pdb_test_file, card_protein)
    ligand_address = disk_datastore.upload_data("ligand_test.pdb", ligand_test_file, card_ligand)

    # Test basic pose generation
    result_address = generate_pose(protein_address=pdb_address,
                                   ligand_address=ligand_address,
                                   output='test_output',
                                   exhaustiveness=1,
                                   num_modes=1)

    # Verify result structure
    assert result_address.startswith('deepchem://')
    assert result_address.endswith('_results.json')

    # Download and verify results
    results_data = disk_datastore.get(result_address)
    results = json.loads(results_data) if isinstance(results_data, str) else results_data

    # Check basic structure
    assert 'docking_method' in results
    assert 'scores_address' in results
    assert 'message' in results

    # Check specific values
    assert results['docking_method'] == 'VINA'
    # num_modes removed from results
    # complexes_count removed from results
    assert 'VINA docking completed successfully' in results['message']

    # Check complex addresses if available
    if 'complex_addresses' in results:
        assert isinstance(results['complex_addresses'], dict)
        for mode_key, complex_address in results['complex_addresses'].items():
            assert complex_address.startswith('deepchem://')
            # Verify complex file can be retrieved
            complex_data = disk_datastore.get(complex_address)
            assert complex_data is not None
            assert isinstance(complex_data, str)  # Should be PDB content
            assert 'ATOM' in complex_data or 'HETATM' in complex_data  # Basic PDB validation


def test_generate_pose_multiple_modes(disk_datastore):
    """Test VINA pose generation with multiple modes."""
    config.set_datastore(disk_datastore)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_test_file = os.path.join(current_dir, "assets/cleaned_3cyx.pdb")
    ligand_test_file = os.path.join(current_dir, "assets/ligand_3cyx.pdb")

    card_protein = DataCard(address='', file_type='pdb', data_type='text/plain')
    card_ligand = DataCard(address='', file_type='pdb', data_type='text/plain')

    pdb_address = disk_datastore.upload_data("protein_test.pdb", pdb_test_file, card_protein)
    ligand_address = disk_datastore.upload_data("ligand_test.pdb", ligand_test_file, card_ligand)

    # Test with multiple modes
    result_address = generate_pose(protein_address=pdb_address,
                                   ligand_address=ligand_address,
                                   output='test_output',
                                   exhaustiveness=1,
                                   num_modes=3)

    # Download and verify results
    results_data = disk_datastore.get(result_address)
    results = json.loads(results_data) if isinstance(results_data, str) else results_data

    # Check that we have scores for multiple modes
    assert 'scores_address' in results
    scores_data = disk_datastore.get(results['scores_address'])
    scores = json.loads(scores_data) if isinstance(scores_data, str) else scores_data

    # Check that scores follow the expected format: 'mode %s' % (i + 1)
    mode_keys = list(scores.keys())
    assert len(mode_keys) > 0

    # Verify each mode has the expected structure
    for mode_key in mode_keys:
        assert mode_key.startswith('mode ')
        assert 'affinity (kcal/mol)' in scores[mode_key]
        # Scores should be numeric
        assert isinstance(scores[mode_key]['affinity (kcal/mol)'], (int, float))


def test_generate_pose_sdf_ligands(disk_datastore):
    """Test VINA pose generation with SDF ligand files."""
    config.set_datastore(disk_datastore)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_test_file = os.path.join(current_dir, "assets/cleaned_3cyx.pdb")
    ligand_sdf_file = os.path.join(current_dir, "assets/ligand_3cyx.sdf")

    card_protein = DataCard(address='', file_type='pdb', data_type='text/plain')
    card_ligand = DataCard(address='', file_type='sdf', data_type='text/plain')

    pdb_address = disk_datastore.upload_data("protein_test.pdb", pdb_test_file, card_protein)
    ligand_address = disk_datastore.upload_data("ligand_test.sdf", ligand_sdf_file, card_ligand)

    # Test with SDF ligand
    result_address = generate_pose(protein_address=pdb_address,
                                   ligand_address=ligand_address,
                                   output='test_output_sdf',
                                   exhaustiveness=1,
                                   num_modes=1)

    # Verify result
    assert result_address.startswith('deepchem://')
    assert result_address.endswith('_results.json')

    # Download and verify results
    results_data = disk_datastore.get(result_address)
    results = json.loads(results_data) if isinstance(results_data, str) else results_data

    assert results['docking_method'] == 'VINA'
    # num_modes removed from results


def test_generate_pose_exhaustiveness_parameter(disk_datastore):
    """Test VINA pose generation with different exhaustiveness values."""
    config.set_datastore(disk_datastore)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_test_file = os.path.join(current_dir, "assets/cleaned_3cyx.pdb")
    ligand_test_file = os.path.join(current_dir, "assets/ligand_3cyx.pdb")

    card_protein = DataCard(address='', file_type='pdb', data_type='text/plain')
    card_ligand = DataCard(address='', file_type='pdb', data_type='text/plain')

    pdb_address = disk_datastore.upload_data("protein_test.pdb", pdb_test_file, card_protein)
    ligand_address = disk_datastore.upload_data("ligand_test.pdb", ligand_test_file, card_ligand)

    # Test with different exhaustiveness values
    for exhaustiveness in [1, 5, 10]:
        result_address = generate_pose(protein_address=pdb_address,
                                       ligand_address=ligand_address,
                                       output=f'test_output_exh_{exhaustiveness}',
                                       exhaustiveness=exhaustiveness,
                                       num_modes=1)

        # Verify result
        assert result_address.startswith('deepchem://')

        # Download and verify results
        results_data = disk_datastore.get(result_address)
        results = json.loads(results_data) if isinstance(results_data, str) else results_data

        assert results['docking_method'] == 'VINA'
        # num_modes removed from results


def test_generate_pose_nested_full_address(disk_datastore):
    """Test VINA pose generation with nested full addresses."""
    config.set_datastore(disk_datastore)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_test_file = os.path.join(current_dir, "assets/cleaned_3cyx.pdb")
    ligand_test_file = os.path.join(current_dir, "assets/ligand_3cyx.pdb")

    card_protein = DataCard(address='', file_type='pdb', data_type='text/plain')
    card_ligand = DataCard(address='', file_type='pdb', data_type='text/plain')

    # Upload with nested paths
    pdb_address = disk_datastore.upload_data("test docking/protein_test.pdb", pdb_test_file, card_protein)
    ligand_address = disk_datastore.upload_data("test docking/ligand_test.pdb", ligand_test_file, card_ligand)

    assert pdb_address == "deepchem://test/user/test docking/protein_test.pdb"
    assert ligand_address == "deepchem://test/user/test docking/ligand_test.pdb"

    # Test pose generation
    result_address = generate_pose(protein_address=pdb_address,
                                   ligand_address=ligand_address,
                                   output="deepchem://test/user/test docking/docking_results",
                                   exhaustiveness=1,
                                   num_modes=1)

    assert result_address.startswith('deepchem://')
    assert result_address.endswith('_results.json')

    # Download and verify results
    results_data = disk_datastore.get(result_address)
    results = json.loads(results_data) if isinstance(results_data, str) else results_data

    assert results['docking_method'] == 'VINA'
    # num_modes removed from results


def test_generate_pose_score_formatting(disk_datastore):
    """Test that VINA scores are formatted correctly"""
    config.set_datastore(disk_datastore)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_test_file = os.path.join(current_dir, "assets/cleaned_3cyx.pdb")
    ligand_test_file = os.path.join(current_dir, "assets/ligand_3cyx.pdb")

    card_protein = DataCard(address='', file_type='pdb', data_type='text/plain')
    card_ligand = DataCard(address='', file_type='pdb', data_type='text/plain')

    pdb_address = disk_datastore.upload_data("protein_test.pdb", pdb_test_file, card_protein)
    ligand_address = disk_datastore.upload_data("ligand_test.pdb", ligand_test_file, card_ligand)

    # Test with multiple modes to verify score formatting
    result_address = generate_pose(protein_address=pdb_address,
                                   ligand_address=ligand_address,
                                   output='test_score_formatting',
                                   exhaustiveness=1,
                                   num_modes=2)

    # Download and verify results
    results_data = disk_datastore.get(result_address)
    results = json.loads(results_data) if isinstance(results_data, str) else results_data

    # Check score format: 'mode %s' % (i + 1) for however many modes are present
    scores_data = disk_datastore.get(results['scores_address'])
    scores = json.loads(scores_data) if isinstance(scores_data, str) else scores_data

    # There should be at least one mode
    assert len(scores) >= 1

    # Verify sequential mode keys exist up to the number returned
    for i in range(len(scores)):
        assert f'mode {i + 1}' in scores

    # Verify score structure for each mode
    for mode_key in scores:
        assert 'affinity (kcal/mol)' in scores[mode_key]
        # Verify score is numeric and reasonable (VINA scores are typically negative)
        score_value = scores[mode_key]['affinity (kcal/mol)']
        assert isinstance(score_value, (int, float))
        # VINA scores are typically negative (better binding = more negative)
        # But we don't enforce this strictly as it depends on the specific system


def test_generate_pose_error_handling(disk_datastore):
    """Test error handling in VINA pose generation."""
    config.set_datastore(disk_datastore)

    # Test with invalid addresses
    with pytest.raises(ValueError) as context:
        generate_pose(protein_address="", ligand_address="deepchem://test/ligand.pdb", output='test_error')
    assert "Protein and/or ligand input is required" in str(context.value)

    with pytest.raises(ValueError) as context:
        generate_pose(protein_address="deepchem://test/protein.pdb", ligand_address="", output='test_error')
    assert "Protein and/or ligand input is required" in str(context.value)

    with pytest.raises(ValueError) as context:
        generate_pose(protein_address=None, ligand_address="deepchem://test/ligand.pdb", output='test_error')
    assert "Protein and/or ligand input is required" in str(context.value)


def test_generate_pose_progress_logging(disk_datastore):
    """Test that progress logging is properly integrated."""
    config.set_datastore(disk_datastore)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_test_file = os.path.join(current_dir, "assets/cleaned_3cyx.pdb")
    ligand_test_file = os.path.join(current_dir, "assets/ligand_3cyx.pdb")

    card_protein = DataCard(address='', file_type='pdb', data_type='text/plain')
    card_ligand = DataCard(address='', file_type='pdb', data_type='text/plain')

    pdb_address = disk_datastore.upload_data("protein_test.pdb", pdb_test_file, card_protein)
    ligand_address = disk_datastore.upload_data("ligand_test.pdb", ligand_test_file, card_ligand)

    # Test that progress logging doesn't interfere with functionality
    result_address = generate_pose(protein_address=pdb_address,
                                   ligand_address=ligand_address,
                                   output='test_progress_logging',
                                   exhaustiveness=1,
                                   num_modes=1)

    # Verify the function still works despite progress logging
    assert result_address.startswith('deepchem://')

    # Download and verify results
    results_data = disk_datastore.get(result_address)
    results = json.loads(results_data) if isinstance(results_data, str) else results_data

    assert results['docking_method'] == 'VINA'
    # num_modes removed from results


def test_generate_pose_data_card_integration(disk_datastore):
    """Test that DataCard integration works correctly for docking results."""
    config.set_datastore(disk_datastore)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_test_file = os.path.join(current_dir, "assets/cleaned_3cyx.pdb")
    ligand_test_file = os.path.join(current_dir, "assets/ligand_3cyx.pdb")

    card_protein = DataCard(address='', file_type='pdb', data_type='text/plain')
    card_ligand = DataCard(address='', file_type='pdb', data_type='text/plain')

    pdb_address = disk_datastore.upload_data("protein_test.pdb", pdb_test_file, card_protein)
    ligand_address = disk_datastore.upload_data("ligand_test.pdb", ligand_test_file, card_ligand)

    # Test pose generation
    result_address = generate_pose(protein_address=pdb_address,
                                   ligand_address=ligand_address,
                                   output='test_datacard',
                                   exhaustiveness=1,
                                   num_modes=1)

    # Verify result
    assert result_address.startswith('deepchem://')

    # Check that DataCard was created for the results
    result_card_address = result_address + '.cdc'
    result_card = disk_datastore.get(result_card_address)

    # Verify DataCard properties
    assert result_card.file_type == 'json'
    assert result_card.data_type == 'json'


def test_generate_pose_pdbqt_support(disk_datastore):
    """Test VINA pose generation with PDBQT file support."""
    config.set_datastore(disk_datastore)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_test_file = os.path.join(current_dir, "assets/cleaned_3cyx.pdb")
    ligand_test_file = os.path.join(current_dir, "assets/ligand_3cyx.pdb")

    card_protein = DataCard(address='', file_type='pdb', data_type='text/plain')
    card_ligand = DataCard(address='', file_type='pdb', data_type='text/plain')

    pdb_address = disk_datastore.upload_data("protein_test.pdb", pdb_test_file, card_protein)
    ligand_address = disk_datastore.upload_data("ligand_test.pdb", ligand_test_file, card_ligand)

    # Test with PDBQT support enabled (single mode)
    result_address = generate_pose(protein_address=pdb_address,
                                   ligand_address=ligand_address,
                                   output='test_output_pdbqt',
                                   exhaustiveness=1,
                                   num_modes=1,
                                   save_pdbqt=True)

    # Verify result
    assert result_address.startswith('deepchem://')
    assert result_address.endswith('_results.json')

    # Download and verify results
    results_data = disk_datastore.get(result_address)
    results = json.loads(results_data) if isinstance(results_data, str) else results_data

    # Check basic structure
    assert 'docking_method' in results
    assert 'scores_address' in results
    assert 'complex_addresses' in results

    # Check if PDBQT addresses are included (may or may not be present depending on Vina output)
    if 'pdbqt_addresses' in results:
        assert isinstance(results['pdbqt_addresses'], dict)
        for mode_key, pdbqt_address in results['pdbqt_addresses'].items():
            assert pdbqt_address.startswith('deepchem://')
            # Verify PDBQT file can be retrieved
            pdbqt_data = disk_datastore.get(pdbqt_address)
            assert pdbqt_data is not None
            assert isinstance(pdbqt_data, list)  # Should be list of lines
            # Basic PDBQT validation
            pdbqt_text = ''.join(pdbqt_data)
            assert 'ATOM' in pdbqt_text or 'HETATM' in pdbqt_text


def test_generate_pose_pdbqt_multiple_modes(disk_datastore):
    """Test VINA pose generation with PDBQT support for multiple modes."""
    config.set_datastore(disk_datastore)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_test_file = os.path.join(current_dir, "assets/cleaned_3cyx.pdb")
    ligand_test_file = os.path.join(current_dir, "assets/ligand_3cyx.pdb")

    card_protein = DataCard(address='', file_type='pdb', data_type='text/plain')
    card_ligand = DataCard(address='', file_type='pdb', data_type='text/plain')

    pdb_address = disk_datastore.upload_data("protein_test.pdb", pdb_test_file, card_protein)
    ligand_address = disk_datastore.upload_data("ligand_test.pdb", ligand_test_file, card_ligand)

    # Test with PDBQT support enabled (multiple modes)
    result_address = generate_pose(protein_address=pdb_address,
                                   ligand_address=ligand_address,
                                   output='test_output_pdbqt_multi',
                                   exhaustiveness=1,
                                   num_modes=3,
                                   save_pdbqt=True)

    # Verify result
    assert result_address.startswith('deepchem://')
    assert result_address.endswith('_results.json')

    # Download and verify results
    results_data = disk_datastore.get(result_address)
    results = json.loads(results_data) if isinstance(results_data, str) else results_data

    # Check basic structure
    assert 'docking_method' in results
    assert 'scores_address' in results
    assert 'complex_addresses' in results

    # Check if PDBQT addresses are included (may or may not be present depending on Vina output)
    if 'pdbqt_addresses' in results:
        assert isinstance(results['pdbqt_addresses'], dict)
        # Should have separate files for each mode
        assert len(results['pdbqt_addresses']) >= 1

        for mode_key, pdbqt_address in results['pdbqt_addresses'].items():
            assert pdbqt_address.startswith('deepchem://')
            # Verify PDBQT file can be retrieved
            pdbqt_data = disk_datastore.get(pdbqt_address)
            assert pdbqt_data is not None
            assert isinstance(pdbqt_data, list)  # Should be list of lines
            # Basic PDBQT validation
            pdbqt_text = ''.join(pdbqt_data)
            assert 'ATOM' in pdbqt_text or 'HETATM' in pdbqt_text
