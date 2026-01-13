"""
Unit tests for RelativeBindingFreeEnergy primitive.
"""
from pathlib import Path
from pyds.data import Data
from pyds.primitives.rbfe import RelativeBindingFreeEnergy
from pyds.settings import Settings


ASSETS_DIR = Path(__file__).resolve().parent / "assets"


class TestRelativeBindingFreeEnergy:
    """Unit tests for RelativeBindingFreeEnergy primitive."""

    def test_init(self, test_settings: Settings) -> None:
        """Test RelativeBindingFreeEnergy initialization."""
        client = RelativeBindingFreeEnergy(settings=test_settings)

        assert client.settings == test_settings
        assert client.base_url == "http://localhost:8000"

    def test_run_with_all_parameters(self, test_settings: Settings, live_data_client: Data) -> None:
        """Test RBFE run with all optional parameters."""

        client = RelativeBindingFreeEnergy(settings=test_settings)

        solvent = {"positive_ion": "Na+", "negative_ion": "Cl-", "ion_concentration": 0.1, "neutralize": True}
        rbfe_settings = {"equilibration_length": "2 ns", "production_length": "20 ns"}

        upload_ligands_sdf_address = live_data_client.upload_data(
            file_path=ASSETS_DIR / "1jld_ligand.sdf",
            filename="test_ligands.sdf",
            description="Test ligands SDF file for relative binding free energy",
        )
        assert "dataset_address" in upload_ligands_sdf_address

        upload_cleaned_protein_pdb_address = live_data_client.upload_data(
            file_path=ASSETS_DIR / "1jld_protein.pdb",
            filename="test_protein.pdb",
            description="Test protein PDB file for relative binding free energy",
        )
        assert "dataset_address" in upload_cleaned_protein_pdb_address

        result = client.run(ligands_sdf_address=upload_ligands_sdf_address["dataset_address"],
                            cleaned_protein_pdb_address=upload_cleaned_protein_pdb_address["dataset_address"],
                            solvent=solvent,
                            overridden_rbfe_settings=rbfe_settings,
                            radial_network_central_ligand=None,
                            dry_run=True,
                            run_edges_in_parallel=True,
                            output_key="custom_output",
                            network_type="RADIAL",
                            scorer_type="LOMAP")

        assert 'relative_binding_free_energy_results_address' in result
