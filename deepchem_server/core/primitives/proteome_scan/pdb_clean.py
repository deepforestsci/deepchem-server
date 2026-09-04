"""
pdb_clean primitive.

Given a gene name and UniProt entry ID, this primitive:
1. Fetches all experimental PDB structures from PDBe/EMBL-EBI.
2. Selects the optimal non-overlapping set (best resolution + coverage).
3. Downloads each selected PDB from PDBe.
4. Cleans each PDB with PDBFixer/OpenMM (remove heterogens/water, add H at pH 7).
5. Uploads all cleaned PDBs and a metadata CSV to the DeepChem datastore.
6. Returns the datastore address of a JSON summary.
"""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

from deepchem_server.core.common import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.common.progress_logger import log_progress
from deepchem_server.core.primitives.proteome_scan import cache


class _Range:

    def __init__(self, id: str, start: int, end: int, error_score: float) -> None:
        self.id = id
        self.start = start
        self.end = end
        self.error_score = error_score
        self.coverage = end - start

    def __repr__(self) -> str:
        return f"Range({self.id}, {self.start}, {self.end}, {self.error_score}, coverage={self.coverage})"


def _get_pdbs_df(entry_id: str) -> tuple[pd.DataFrame, int]:
    """Fetch PDB metadata for a UniProt entry from PDBe.

    Parameters
    ----------
    entry_id : str
        UniProt accession (e.g. 'P02144').

    Returns
    -------
    tuple[pd.DataFrame, int]
        DataFrame with columns id, resolution, chains, chain_type, chain_start,
        chain_end, coverage (index set to id), and the sequence length of the
        UniProt entry.
    """
    resp = requests.get(
        f"https://www.ebi.ac.uk/pdbe/api/v2/uniprot/unipdb/{entry_id}",
        headers={"Accept": "application/json"},
        timeout=60,
    )
    resp.raise_for_status()
    data = resp.json()
    entry_data = data[entry_id]
    seq_length: int = entry_data["length"]

    def _get_pdb(entry: dict) -> Optional[dict]:
        if not isinstance(entry, dict):
            return None
        output: dict = {"type": "PDB"}

        pdb_id = str(entry["name"]).upper()
        output["id"] = pdb_id

        output["resolution"] = entry.get("additionalData", {}).get("resolution") or None

        # Compute the mapping range from non-mutation residue entries only.
        # Mutation entries have a "mutation" key; skip them so they don't
        # distort start/end.
        main_residues = [r for r in entry.get("residues", []) if not r.get("mutation")]
        if not main_residues:
            return None
        start = min(r["startIndex"] for r in main_residues)
        end = max(r["endIndex"] for r in main_residues)

        # bestChainId is a first-class field in v2 — no regex needed.
        best_chain = entry.get("bestChainId")
        if not best_chain:
            return None

        # Fetch all chains belonging to this protein entity so we keep every
        # chain produced by the gene, not just the best one.
        try:
            mol_resp = requests.get(
                f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/molecules/{pdb_id.lower()}",
                headers={"Accept": "application/json"},
                timeout=60,
            )
            mol_data = mol_resp.json()
            pdb_key = pdb_id.lower()
            all_entity_chains = [best_chain]
            for mol in mol_data.get(pdb_key, []):
                if mol.get("molecule_type") == "polypeptide(L)" and best_chain in mol.get("in_chains", []):
                    all_entity_chains = mol["in_chains"]
                    break
            output["chains"] = "/".join(all_entity_chains) + f"={start}-{end}"
        except Exception:
            output["chains"] = best_chain + f"={start}-{end}"

        return output

    pdb_list_raw = entry_data["data"]
    with ThreadPoolExecutor(max_workers=5) as executor:
        pdb_list = list(executor.map(_get_pdb, pdb_list_raw))
    pdb_list = [p for p in pdb_list if isinstance(p, dict)]

    pdbs_df = pd.DataFrame(pdb_list)
    pdbs_df["id"] = pdbs_df["id"].str.upper()
    pdbs_df[["chain_type", "chain_start", "chain_end"]] = pdbs_df["chains"].str.extract(r"(.+)?=(\d+)-(\d+)")
    pdbs_df = pdbs_df.dropna(subset=["chain_type", "chain_start", "chain_end"])
    pdbs_df["chain_type"] = pdbs_df["chain_type"].apply(lambda x: str(x).split("/"))
    pdbs_df["chain_start"] = pd.to_numeric(pdbs_df["chain_start"], errors="coerce")
    pdbs_df["chain_end"] = pd.to_numeric(pdbs_df["chain_end"], errors="coerce")
    pdbs_df = pdbs_df.dropna(subset=["chain_start", "chain_end"])
    pdbs_df["chain_start"] = pdbs_df["chain_start"].astype(int)
    pdbs_df["chain_end"] = pdbs_df["chain_end"].astype(int)
    pdbs_df["coverage"] = pdbs_df["chain_end"] - pdbs_df["chain_start"]
    pdbs_df = pdbs_df.set_index("id", drop=False)
    return pdbs_df, seq_length


def _select_ranges(ranges: list[_Range]) -> list[_Range]:
    """Greedy selection of a non-overlapping set of sequence ranges.

    Sorts by (start, end) and then iterates. Non-overlapping ranges are always
    kept. For overlapping ones, the current range replaces the previous
    selection when it has better resolution or similar resolution with
    greater coverage.

    Parameters
    ----------
    ranges : list[_Range]
        Input range objects.

    Returns
    -------
    list[_Range]
        Selected non-overlapping ranges.
    """
    sorted_ranges = sorted(ranges, key=lambda r: (r.start, r.end))
    selected: list[_Range] = []
    last_end = float("-inf")

    for cur in sorted_ranges:
        if cur.start > last_end:
            selected.append(cur)
            last_end = cur.end
        else:
            if (cur.error_score < selected[-1].error_score or abs(cur.error_score - selected[-1].error_score) < 0.5):
                if cur.coverage > selected[-1].coverage:
                    selected.append(cur)
                last_end = max(last_end, cur.end)

    return selected


def _get_optimal_pdbs_df(df: pd.DataFrame, seq_length: int, min_res_val: float = 2.5) -> pd.DataFrame:
    """Select the best PDB structures for a given protein.

    Builds two candidate sets — one filtered by the resolution threshold and
    one by best coverage — runs greedy range selection on each, and returns
    the union.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame of PDB candidates with resolution and coverage data.
    seq_length : int
        Full sequence length of the protein.
    min_res_val : float, optional
        Maximum resolution (Angstrom) for the high-res candidate set.
        Default is 2.5.

    Returns
    -------
    pd.DataFrame
        Subset of df containing only the selected optimal structures.
    """
    df = df.dropna(subset=["resolution"])
    df["resolution"] = pd.to_numeric(df["resolution"], errors="coerce")
    df = df.dropna(subset=["resolution"])

    df["range_obj"] = df.apply(lambda x: _Range(x.id, x.chain_start, x.chain_end, x.resolution), axis=1)
    min_res_for_max_cov = float(min(df[df["coverage"] == max(df["coverage"])]["resolution"].values))

    df_1 = df[df["resolution"] <= min_res_val]
    df_2 = df[df["resolution"] <= min_res_for_max_cov]

    selected_1 = _select_ranges(df_1["range_obj"].to_list())
    selected_2 = _select_ranges(df_2["range_obj"].to_list())

    selected = list(set(selected_1 + selected_2))
    ids = [r.id for r in selected]
    return df.loc[ids]


def _download_pdb(gene_name: str, pdb_id: str, scan_id: str) -> bool:
    """Download a PDB file from PDBe into the local cache.

    Parameters
    ----------
    gene_name : str
        Gene name (used to construct the cache path).
    pdb_id : str
        PDB identifier (uppercase).
    scan_id : str
        Scan run identifier.

    Returns
    -------
    bool
        True if the download succeeded, False otherwise.
    """
    out_path = cache.pdb_raw_path(scan_id, gene_name, pdb_id)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists():
        return True

    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pdb_id.lower()}.ent"
    try:
        resp = requests.get(url, timeout=60)
        if resp.status_code != 200 or not resp.content:
            raise RuntimeError(f"bad response status {resp.status_code}")
        out_path.write_bytes(resp.content)
        return True
    except Exception:
        if out_path.exists():
            out_path.unlink()
        return False


def _get_chain_ids(pdb_path: Path) -> list[str]:
    """Return all chain IDs present in a PDB file.

    Parameters
    ----------
    pdb_path : Path
        Path to the PDB file.

    Returns
    -------
    list[str]
        Chain IDs found in the structure.
    """
    from Bio import PDB  # type: ignore

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", str(pdb_path))
    chain_ids: set[str] = set()
    for model in structure:
        for chain in model:
            chain_ids.add(chain.id)
    return list(chain_ids)


def _pdb_cleaner(gene_name: str,
                 pdb_id: str,
                 scan_id: str,
                 remove_chains: list[str],
                 keep_water: bool = False,
                 add_missing_hydrogens_pH: float = 7.0) -> Optional[Path]:
    """Clean a downloaded PDB with PDBFixer and OpenMM.

    Removes specified chains, strips heterogens and water, then adds missing
    hydrogens at pH 7.0.

    Parameters
    ----------
    gene_name : str
        Gene name.
    pdb_id : str
        PDB identifier (uppercase).
    scan_id : str
        Scan run identifier.
    remove_chains : list[str]
        Chain IDs to remove before cleaning.
    keep_water : bool
        Whether to keep water molecules while removing heterogens. True by default.
    add_missing_hydrogens_pH : float
        The pH based on which to select hydrogens. 7.0 by default.

    Returns
    -------
    Path or None
        Path to the cleaned PDB file, or None if cleaning failed.

    Raises
    ------
    ImportError
        If pdbfixer or openmm is not installed.
    """
    try:
        from pdbfixer.pdbfixer import PDBFixer  # type: ignore
    except ImportError as exc:
        raise ImportError("pdbfixer is required for pdb_clean but not installed") from exc
    try:
        from openmm.app import PDBFile  # type: ignore
    except ImportError as exc:
        raise ImportError("openmm is required for pdb_clean but not installed") from exc

    pdb_path = cache.pdb_raw_path(scan_id, gene_name, pdb_id)
    output_path = cache.pdb_clean_path(scan_id, gene_name, pdb_id)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        fixer = PDBFixer(str(pdb_path))
        if remove_chains:
            fixer.removeChains(chainIds=remove_chains)
        fixer.removeHeterogens(keepWater=keep_water)
        try:
            fixer.addMissingHydrogens(pH=add_missing_hydrogens_pH)
        except Exception:
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(pH=add_missing_hydrogens_pH)
        with open(output_path, "w") as fp:
            PDBFile.writeFile(fixer.topology, fixer.positions, fp)
        return output_path
    except Exception as exc:
        log_progress("pdb_clean", 75, f"pdbfixer failed for {pdb_id}: {exc}")
        return None


def pdb_clean(
    gene_name: str,
    entry_id: str,
    scan_id: str,
    output: str,
    min_res_val: float = 2.5,
) -> str:
    """Download and clean the best experimental PDB structures for a gene.

    Fetches all experimental structures for the given UniProt entry, selects
    the optimal non-overlapping coverage set by resolution and sequence
    coverage, cleans each with PDBFixer/OpenMM, uploads everything to the
    datastore, and returns the address of a JSON summary.

    Parameters
    ----------
    gene_name : str
        Human gene name (e.g. 'MAP2K1'). Used to name cache artifacts.
    entry_id : str
        UniProt accession (e.g. 'P02144').
    scan_id : str
        Identifier for the current scan run. All cache artifacts are
        namespaced under this ID.
    output : str
        Prefix used when naming files in the datastore.
    min_res_val : float, optional
        Resolution threshold (Angstrom) for the high-res candidate set.
        Default is 2.5.

    Returns
    -------
    str
        Datastore address of a JSON summary containing gene_name, entry_id,
        scan_id, min_res_val, selected_pdb_ids, cleaned_pdb_addresses, and
        pdbs_metadata_csv_address.

    Raises
    ------
    ValueError
        If the datastore is not configured.
    OverflowError
        If all candidate PDB structures fail to download or clean.
    RuntimeError
        If no cleaned PDB files could be produced.

    Examples
    --------
    >>> addr = pdb_clean("MB", "P02144", "scan_001", "mb_output")
    >>> print(addr)
    deepchem://profile/project/mb_output/MB_pdb_clean_summary.json
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not configured; call config.set_datastore() first.")

    log_progress("pdb_clean", 10, f"fetching PDB metadata for entry {entry_id}")
    pdbs_df, seq_length = _get_pdbs_df(entry_id)

    log_progress("pdb_clean", 30, f"selecting optimal PDB structures (min_res_val={min_res_val})")

    failed_pdbs: Optional[list[str]] = None
    selected_pdbs: pd.DataFrame = pd.DataFrame()
    id_chain_details: dict[str, tuple[list[str], list[str]]] = {}

    while failed_pdbs != []:
        if failed_pdbs is not None:
            pdbs_df = pdbs_df.drop(index=failed_pdbs)
        if len(pdbs_df) == 0:
            raise OverflowError("no pdbs left to check, likely over sized pdbs")

        selected_pdbs = _get_optimal_pdbs_df(pdbs_df, seq_length, min_res_val)
        id_chain_map: dict[str, list[str]] = selected_pdbs["chain_type"].to_dict()

        log_progress("pdb_clean", 50, f"downloading {len(id_chain_map)} PDB files")
        failed_pdbs = [pdb_id for pdb_id in id_chain_map if not _download_pdb(gene_name, pdb_id, scan_id)]
        if failed_pdbs:
            log_progress("pdb_clean", 55, f"failed downloads: {failed_pdbs}, retrying selection")
            continue

        id_chain_details = {}
        for pdb_id, keep_chains in id_chain_map.items():
            raw_path = cache.pdb_raw_path(scan_id, gene_name, pdb_id)
            try:
                all_chains = _get_chain_ids(raw_path)
            except Exception:
                all_chains = keep_chains
            chains_to_remove = list(set(all_chains) - set(keep_chains))
            id_chain_details[pdb_id] = (keep_chains, chains_to_remove)

        log_progress("pdb_clean", 70, "cleaning PDB files with PDBFixer")
        uncleanable: list[str] = []
        for pdb_id, (_, remove_chains) in id_chain_details.items():
            if _pdb_cleaner(gene_name, pdb_id, scan_id, remove_chains) is None:
                uncleanable.append(pdb_id)
        failed_pdbs = uncleanable
        if failed_pdbs:
            log_progress("pdb_clean", 72, f"uncleanable PDBs: {failed_pdbs}, retrying selection")

    log_progress("pdb_clean", 90, "uploading cleaned PDBs and metadata to datastore")

    cleaned_pdb_addresses: dict[str, str] = {}
    for pdb_id in id_chain_details:
        clean_path = cache.pdb_clean_path(scan_id, gene_name, pdb_id)
        if not clean_path.exists():
            continue
        ds_filename = f"{output}/cleaned_g_{gene_name}_p_{pdb_id}.pdb"
        card = DataCard(address="", file_type="pdb", data_type="text/plain")
        addr = datastore.upload_data(
            datastore_filename=ds_filename,
            filename=str(clean_path),
            card=card,
        )
        if addr:
            cleaned_pdb_addresses[pdb_id] = addr

    if not cleaned_pdb_addresses:
        raise RuntimeError(f"No cleaned PDB files could be produced for gene {gene_name} (entry {entry_id})")

    csv_filename = f"{output}/{gene_name}_pdbs.csv"
    csv_card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    csv_df = selected_pdbs.drop(columns=["range_obj"], errors="ignore")
    pdbs_metadata_csv_address = datastore.upload_data_from_memory(
        data=csv_df,
        datastore_filename=csv_filename,
        card=csv_card,
    )

    summary = {
        "gene_name": gene_name,
        "entry_id": entry_id,
        "scan_id": scan_id,
        "min_res_val": min_res_val,
        "selected_pdb_ids": list(cleaned_pdb_addresses.keys()),
        "cleaned_pdb_addresses": cleaned_pdb_addresses,
        "pdbs_metadata_csv_address": pdbs_metadata_csv_address,
    }
    summary_filename = f"{output}/{gene_name}_pdb_clean_summary.json"
    summary_card = DataCard(address="", file_type="json", data_type="json")
    result_address = datastore.upload_data_from_memory(
        data=json.dumps(summary),
        datastore_filename=summary_filename,
        card=summary_card,
    )

    if result_address is None:
        raise ValueError("Failed to upload pdb_clean results summary to datastore")

    log_progress(
        "pdb_clean",
        100,
        f"pdb_clean completed: {len(cleaned_pdb_addresses)} cleaned PDBs for {gene_name}",
    )
    return result_address
