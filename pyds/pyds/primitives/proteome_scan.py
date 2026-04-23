"""
ProteomeScan primitives module for DeepChem Server (client side).

Contains client wrappers for the proteome scan pipeline primitives:

- PdbClean: fetch and clean PDBs for a gene
- RunDocking: run VINA docking for a (gene, ligand) pair
- ParseResults: aggregate per-ligand docking results
- RunMultiPoseAnalysis: pose-binding analysis for docked complexes

The module also exposes the high-level ProteomeScan pipeline
client, which orchestrates the four primitives end-to-end.
"""

from typing import Any, Dict, List, Optional

from pyds.settings import Settings

from .base import Primitive


class PdbClean(Primitive):
    """Client for /primitive/proteome-scan/pdb-clean."""

    def run(
        self,
        gene_name: str,
        entry_id: str,
        scan_id: str,
        output: str,
        min_res_val: float = 2.5,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        profile, project = self.validate_common_params(profile_name, project_name)

        payload = {
            "profile_name": profile,
            "project_name": project,
            "gene_name": gene_name,
            "entry_id": entry_id,
            "scan_id": scan_id,
            "output": output,
            "min_res_val": min_res_val,
        }
        response = self._post("/primitive/proteome-scan/pdb-clean", json=payload)
        return self._validate_response(response)


class FilterPromiscuousTargets(Primitive):
    """
    Primitive for filtering promiscuous targets from docking scan results.

    This class submits jobs to the DeepChem Server API's
    /primitive/filter-promiscuous-targets endpoint.
    """

    def run(
        self,
        scan_result_addresses: List[str],
        thresholds: List[List[int]],
        output: str,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Run the filter_promiscuous_targets primitive.

        Parameters
        ----------
        scan_result_addresses: List[str]
            DeepChem addresses of per-ligand scan result
                CSV files. Each CSV must contain a gene_name column and rows
                should be sorted from best (top) to worst docking score.
        thresholds: List[List[int]]
            List of [m, n] pairs where m is the percentile cutoff
                (0-100) and n is the minimum occurrence count across ligands.
                For example, [[15, 1]] marks a gene as promiscuous if it
                appears in the top 15% of at least 1 ligand's results.
        output: str
            Output name prefix used for all uploaded result files.
        profile_name: Optional[str]
            Profile name (uses settings if not provided).
        project_name: Optional[str]
            Project name (uses settings if not provided).

        Returns
        -------
            Response dict with filter_results_address pointing to a
            JSON file that contains promiscuous target information.
        """
        profile, project = self.validate_common_params(profile_name, project_name)

        data = {
            "profile_name": profile,
            "project_name": project,
            "scan_result_addresses": scan_result_addresses,
            "thresholds": thresholds,
            "output": output,
        }

        response = self._post("/primitive/filter-promiscuous-targets", json=data)
        return self._validate_response(response)


class LigandPrep(Primitive):
    """
    Primitive for preparing ligands (SMILES to 3D SDF) via the DeepChem Server API.

    Submits jobs to the /primitive/ligand-prep endpoint.
    """

    def run(
        self,
        smiles: str,
        output: str,
        ligand_name: str = "",
        random_seed: int = 42,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Run the ligand preparation primitive.

        Parameters
        ----------
        smiles : str
            Input SMILES string.
        output : str
            Output key for the SDF in the datastore.
        ligand_name : str
            Optional molecule name stored in the SDF.
        random_seed : int
            Seed for 3D embedding.
        profile_name : str or None
            Profile name; uses settings when omitted.
        project_name : str or None
            Project name; uses settings when omitted.

        Returns
        -------
        dict
            Response containing the ligand_sdf_address entry.
        """
        profile, project = self.validate_common_params(profile_name, project_name)

        data: Dict[str, Any] = {
            "profile_name": profile,
            "project_name": project,
            "smiles": smiles,
            "output": output,
            "ligand_name": ligand_name,
            "random_seed": random_seed,
        }

        response = self._post("/primitive/proteome-scan/ligand-prep", json=data)
        return self._validate_response(response)


class ParseResults(Primitive):
    """Client for /primitive/proteome-scan-parse-results."""

    def run(
        self,
        scan_id: str,
        ligands: List[str],
        output: str,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        profile, project = self.validate_common_params(profile_name, project_name)

        payload = {
            "profile_name": profile,
            "project_name": project,
            "scan_id": scan_id,
            "ligands": ligands,
            "output": output,
        }
        response = self._post("/primitive/proteome-scan/parse-results", json=payload)
        return self._validate_response(response)


class RunDocking(Primitive):
    """Client for /primitive/proteome-scan-docking."""

    def run(
        self,
        gene_name: str,
        ligand_name: str,
        ligand_address: str,
        scan_id: str,
        output: str,
        exhaustiveness: int = 32,
        num_modes: int = 8,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        profile, project = self.validate_common_params(profile_name, project_name)

        payload = {
            "profile_name": profile,
            "project_name": project,
            "gene_name": gene_name,
            "ligand_name": ligand_name,
            "ligand_address": ligand_address,
            "scan_id": scan_id,
            "output": output,
            "exhaustiveness": exhaustiveness,
            "num_modes": num_modes,
        }
        response = self._post("/primitive/proteome-scan/docking", json=payload)
        return self._validate_response(response)


class RunMultiPoseAnalysis(Primitive):
    """Client for /primitive/proteome-scan-multi-pose-analysis."""

    def run(
        self,
        complex_addresses: List[str],
        output: str,
        num_processes: int = 4,
        is_clean_up: bool = True,
        scan_id: Optional[str] = None,
        profile_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        profile, project = self.validate_common_params(profile_name, project_name)

        payload = {
            "profile_name": profile,
            "project_name": project,
            "complex_addresses": complex_addresses,
            "output": output,
            "num_processes": num_processes,
            "is_clean_up": is_clean_up,
            "scan_id": scan_id,
        }
        response = self._post("/primitive/proteome-scan/multi-pose-analysis", json=payload)
        return self._validate_response(response)


class ProteomeScan:
    """
    High-level client for running proteome-wide docking scans.

    Orchestrates the ProteomeScan pipeline end-to-end using the
    deepchem-server primitives:

    1. pdb_clean for each gene
    2. run_docking for each (gene, ligand) pair
    3. parse_results to aggregate per-ligand top-score tables
    4. optional run_multi_pose_analysis

    Example
    -------
    >>> from pyds import Settings
    >>> from pyds.primitives.proteome_scan import ProteomeScan
    >>> settings = Settings()
    >>> settings.set_profile("default")
    >>> settings.set_project("my_project")
    >>> scanner = ProteomeScan(settings)
    >>> result = scanner.run_pipeline(
    ...     scan_id="scan_20261231",
    ...     ligands=["Trametinib", "Tucatinib"],
    ...     gene_names=["GBA3", "SLC7A11"],
    ...     gene_entry_map={"GBA3": "Q9H3H0", "SLC7A11": "Q9UPY5"},
    ...     ligand_address_map={
    ...         "Trametinib": "deepchem://default/my_project/trametinib.sdf",
    ...         "Tucatinib": "deepchem://default/my_project/tucatinib.sdf",
    ...     },
    ... )
    """

    def __init__(self, settings: Optional[Settings] = None):
        self.settings = settings or Settings()
        self.pdb_clean = PdbClean(self.settings)
        self.run_docking = RunDocking(self.settings)
        self.parse_results = ParseResults(self.settings)
        self.run_multi_pose_analysis = RunMultiPoseAnalysis(self.settings)

    def run_pipeline(
        self,
        scan_id: str,
        ligands: List[str],
        gene_names: List[str],
        gene_entry_map: Dict[str, str],
        ligand_address_map: Dict[str, str],
        exhaustiveness: int = 32,
        num_modes: int = 8,
        output_prefix: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Run the complete ProteomeScan pipeline.

        Parameters
        ----------
        scan_id: str
            Scan identifier used by all primitives to share on-disk
            artifacts on the server.
        ligands: list of str
            Ligand display names (must match keys in
            ligand_address_map).
        gene_names: list of str
            Gene symbols to scan. Each must appear in
            gene_entry_map.
        gene_entry_map: dict
            Mapping gene_name -> uniprot_entry_id.
        ligand_address_map: dict
            Mapping ligand_name -> datastore address of prepared SDF.
        exhaustiveness, num_modes: int
            Passed to run_docking.
        output_prefix: str, optional
            Prefix used for every uploaded artifact. Defaults to scan_id.

        Returns
        -------
        dict
            {"scan_id": ..., "gene_pdbs": {...},
            "docking_results": {...}, "parsed_results": <address>}
        """
        output_prefix = output_prefix or scan_id

        gene_pdbs: Dict[str, str] = {}
        for gene_name in gene_names:
            entry_id = gene_entry_map.get(gene_name)
            if not entry_id:
                print(f"Warning: no entry_id for {gene_name}, skipping")
                continue
            try:
                result = self.pdb_clean.run(
                    gene_name=gene_name,
                    entry_id=entry_id,
                    scan_id=scan_id,
                    output=f"{output_prefix}_{gene_name}_pdbs",
                )
                gene_pdbs[gene_name] = result["results_address"]
                print(f"cleaned PDBs fetched for {gene_name}: {result['results_address']}")
            except Exception as e:  # noqa: BLE001
                print(f"failed to get cleaned PDBs for {gene_name}: {e}")

        if not gene_pdbs:
            raise ValueError("No gene PDBs could be fetched")

        docking_results: Dict[str, Dict[str, str]] = {}
        for gene_name in gene_pdbs:
            docking_results[gene_name] = {}
            for ligand in ligands:
                ligand_addr = ligand_address_map[ligand]
                try:
                    res = self.run_docking.run(
                        gene_name=gene_name,
                        ligand_name=ligand,
                        ligand_address=ligand_addr,
                        scan_id=scan_id,
                        output=f"{output_prefix}_{gene_name}_{ligand}",
                        exhaustiveness=exhaustiveness,
                        num_modes=num_modes,
                    )
                    docking_results[gene_name][ligand] = res["results_address"]
                    print(f"docking results for {gene_name} {ligand}: {res['results_address']}")
                except Exception as e:  # noqa: BLE001
                    print(f"failed to run docking for {gene_name} {ligand}: {e}")

        parsed = self.parse_results.run(
            scan_id=scan_id,
            ligands=ligands,
            output=f"{output_prefix}_parse",
        )

        return {
            "scan_id": scan_id,
            "gene_pdbs": gene_pdbs,
            "docking_results": docking_results,
            "parsed_results": parsed["results_address"],
        }

    def get_gene_pdb(
        self,
        gene_name: str,
        entry_id: str,
        scan_id: str,
        output: Optional[str] = None,
    ) -> str:
        """Fetch cleaned PDBs for a single gene and return the datastore
        address of the metadata JSON."""
        result = self.pdb_clean.run(
            gene_name=gene_name,
            entry_id=entry_id,
            scan_id=scan_id,
            output=output or f"{scan_id}_{gene_name}_pdbs",
        )
        return result["results_address"]

    def dock_pair(
        self,
        gene_name: str,
        ligand_name: str,
        ligand_address: str,
        scan_id: str,
        output: Optional[str] = None,
        exhaustiveness: int = 32,
        num_modes: int = 8,
    ) -> str:
        """Run docking for a single (gene, ligand) pair and return the
        datastore address of the run summary JSON."""
        result = self.run_docking.run(
            gene_name=gene_name,
            ligand_name=ligand_name,
            ligand_address=ligand_address,
            scan_id=scan_id,
            output=output or f"{scan_id}_{gene_name}_{ligand_name}",
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
        )
        return result["results_address"]
