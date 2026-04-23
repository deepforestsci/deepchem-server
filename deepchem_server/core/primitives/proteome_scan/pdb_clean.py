"""
Fetch and clean PDB structures for a target gene.
"""

import json
import os
import re
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Dict, List, Optional, Tuple

import pandas as pd
import requests

from deepchem_server.core.common import config
from deepchem_server.core.common.cards import DataCard
from deepchem_server.core.common.progress_logger import log_progress
from deepchem_server.core.primitives.proteome_scan import cache as ps_cache


def clean_pdb_file(
    input_pdb_path: str,
    output_pdb_path: str,
    remove_chains: Optional[List[str]] = None,
    remove_heterogens: bool = True,
    remove_water: bool = True,
    add_hydrogens: bool = True,
    ph: float = 7.0,
) -> Optional[str]:
    """
    Clean a local PDB file and write the cleaned output to disk.

    Parameters
    ----------
    input_pdb_path : str
        Path to the input PDB file.
    output_pdb_path : str
        Path where the cleaned PDB will be written.
    remove_chains : list[str] or None, default None
        Chain identifiers to remove before writing the cleaned structure.
    remove_heterogens : bool, default True
        If True, remove heterogens from the structure.
    remove_water : bool, default True
        If True and remove_heterogens is True, water molecules are removed too.
    add_hydrogens : bool, default True
        If True, add missing hydrogens with the provided pH.
    ph : float, default 7.0
        pH value used when adding hydrogens.

    Returns
    -------
    str or None
        Output path if cleaning succeeded, otherwise None.

    Examples
    --------
    >>> out = clean_pdb_file("in.pdb", "out.pdb")
    """
    try:
        from pdbfixer.pdbfixer import PDBFixer  # type: ignore
    except ImportError:
        raise ImportError("pdbfixer is required for pdb_clean but not installed")
    try:
        from openmm.app import PDBFile  # type: ignore
    except ImportError:
        raise ImportError("OpenMM is required for pdb_clean but not installed")

    try:
        fixer = PDBFixer(input_pdb_path)
        if remove_chains:
            fixer.removeChains(chainIds=list(remove_chains))
        if remove_heterogens:
            fixer.removeHeterogens(not remove_water)
        if add_hydrogens:
            try:
                fixer.addMissingHydrogens(ph)
            except Exception:
                fixer.findMissingResidues()
                fixer.findMissingAtoms()
                fixer.addMissingAtoms()
                fixer.addMissingHydrogens(ph)

        os.makedirs(os.path.dirname(output_pdb_path) or ".", exist_ok=True)
        with open(output_pdb_path, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        return output_pdb_path
    except Exception as e:  # noqa: BLE001
        log_progress("pdb_clean", 0, f"failed PDB fixing: {input_pdb_path} => {e}")
        return None


def run_on_multiple_threads(fn: Callable, values: list, max_workers: int) -> list:
    """
    Apply a function over a list using a thread pool.

    Parameters
    ----------
    fn : Callable
        Function to apply to each input value.
    values : list
        Input values to map over.
    max_workers : int
        Maximum number of worker threads.

    Returns
    -------
    list
        Results in the same order as the input values.

    Examples
    --------
    >>> run_on_multiple_threads(lambda x: x + 1, [1, 2, 3], max_workers=2)
    [2, 3, 4]
    """
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(fn, values))
    return list(results)


def get_protein_details(protein_id: str) -> Tuple[int, str]:
    """
    Fetch protein sequence details from the EMBL-EBI proteins API.

    Parameters
    ----------
    protein_id : str
        UniProt protein identifier.

    Returns
    -------
    tuple[int, str]
        Tuple of (sequence_length, sequence_string).

    Examples
    --------
    >>> length, seq = get_protein_details("Q9H3H0")
    """
    last_error: Optional[Exception] = None
    for _ in range(5):
        try:
            search_results = requests.get(
                f"https://www.ebi.ac.uk/proteins/api/proteins/{protein_id}",
                headers={"Accept": "application/json"},
                timeout=60,
            )
            search_results.raise_for_status()
            data = search_results.json()
            if not isinstance(data, dict) or "sequence" not in data:
                raise RuntimeError("unexpected proteins API payload shape")
            length, seq = data["sequence"]["length"], data["sequence"]["sequence"]
            return length, seq
        except Exception as e:  # noqa: BLE001
            last_error = e
    raise RuntimeError(f"failed to fetch protein details for {protein_id}: {last_error}")


def get_pdbs_df(protein_id: str) -> pd.DataFrame:
    """
    Fetch PDB metadata for a UniProt protein identifier.

    Parameters
    ----------
    protein_id : str
        UniProt protein identifier.

    Returns
    -------
    pandas.DataFrame
        Table of candidate PDB structures with coverage, chains, and resolution.

    Examples
    --------
    >>> df = get_pdbs_df("Q9H3H0")
    """
    last_error: Optional[Exception] = None
    data: object = {}
    for _ in range(5):
        try:
            search_results = requests.get(
                f"https://www.ebi.ac.uk/pdbe/graph-api/uniprot/protvista/unipdb/{protein_id}",
                headers={"Accept": "application/json"},
                timeout=60,
            )
            search_results.raise_for_status()
            data = search_results.json()
            if isinstance(data, dict) and protein_id in data:
                break
            raise RuntimeError("unexpected PDBe unipdb payload shape")
        except Exception as e:  # noqa: BLE001
            last_error = e
            data = {}
    if not (isinstance(data, dict) and protein_id in data):
        raise RuntimeError(f"failed to fetch PDB metadata for {protein_id}: {last_error}")

    def get_pdb(pdb_struct):
        if not isinstance(pdb_struct, dict):
            return None
        output = {"type": "PDB"}
        # PDBe/ProtVista can return lowercase ids (e.g. "3rgk"), while other
        # PDBe endpoints key response payloads by uppercase ids (e.g. "3RGK").
        # Normalize early to avoid KeyError mismatches across calls.
        pdb_id = str(pdb_struct["label"]["id"]).upper()
        output["id"] = pdb_id
        output["resolution"] = pdb_struct["label"].get("resolution")
        if not output["resolution"]:
            output["resolution"] = None
        start = pdb_struct["locations"][0]["fragments"][0]["start"]
        end = pdb_struct["locations"][0]["fragments"][-1]["end"]
        pattern = r'href="([^"]+)"'
        match = re.search(pattern, pdb_struct["locations"][0]["fragments"][0]["tooltipContent"])

        if match:

            def get_chains(match):
                href_link = match.group(1)
                chain_search = requests.get(
                    f"https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/protvista/chains/{pdb_id.lower()}/{href_link.split('/')[-1]}"
                )
                chain_data = chain_search.json()
                pdb_key = pdb_id
                if pdb_key not in chain_data:
                    for candidate in (pdb_id.lower(), pdb_id.upper()):
                        if candidate in chain_data:
                            pdb_key = candidate
                            break
                if pdb_key not in chain_data and len(chain_data) == 1:
                    pdb_key = next(iter(chain_data))
                chains = set()
                node = chain_data.get(pdb_key)
                if not isinstance(node, dict):
                    raise RuntimeError("unexpected chains payload shape")
                for data_ in node["tracks"][0]["data"]:
                    text = data_["label"]
                    pat = r"\(auth (.*)\)"
                    m = re.search(pat, text)
                    if m:
                        chains.add(m.group(1))
                return "/".join(chains) + "=" + str(start) + "-" + str(end)

            fetched = False
            count = 0
            while not fetched and count < 5:
                try:
                    output["chains"] = get_chains(match)
                    fetched = True
                except json.decoder.JSONDecodeError:
                    count += 1
                except Exception:  # noqa: BLE001
                    # Best-effort fallback: if chain metadata lookup fails, assume a
                    # single-chain structure so the pipeline can proceed.
                    output["chains"] = f"A={start}-{end}"
                    fetched = True
        else:
            output["chains"] = None
        return output

    raw_items = data[protein_id]["tracks"][0]["data"]
    pdb_list = run_on_multiple_threads(get_pdb, raw_items, max_workers=5)
    pdb_list = [p for p in pdb_list if isinstance(p, dict)]

    pdbs_df = pd.DataFrame(pdb_list)
    pdbs_df["id"] = pdbs_df["id"].str.upper()
    pdbs_df[["chain_type", "chain_start", "chain_end"]] = pdbs_df["chains"].str.extract(
        r"(.+)?=(\d+)-(\d+)"
    )
    pdbs_df = pdbs_df.dropna(subset=["chain_type", "chain_start", "chain_end"])
    pdbs_df["chain_type"] = pdbs_df["chain_type"].apply(lambda x: str(x).split("/"))
    pdbs_df["chain_start"] = pd.to_numeric(pdbs_df["chain_start"], errors="coerce")
    pdbs_df["chain_end"] = pd.to_numeric(pdbs_df["chain_end"], errors="coerce")
    pdbs_df = pdbs_df.dropna(subset=["chain_start", "chain_end"])
    pdbs_df["chain_start"] = pdbs_df["chain_start"].astype(int)
    pdbs_df["chain_end"] = pdbs_df["chain_end"].astype(int)
    pdbs_df["coverage"] = pdbs_df["chain_end"] - pdbs_df["chain_start"]
    pdbs_df = pdbs_df.set_index("id", drop=False)
    return pdbs_df


class Range:
    """
    Range of sequence positions covered by a PDB chain.

    Parameters
    ----------
    id : str
        PDB identifier.
    start : int
        Start position (inclusive).
    end : int
        End position (inclusive).
    error_score : float
        Score used to compare ranges, typically resolution.

    Examples
    --------
    >>> r = Range("1ABC", 1, 100, 2.0)
    >>> r.coverage
    99
    """

    def __init__(self, id: str, start: int, end: int, error_score: float) -> None:
        self.id = id
        self.start = start
        self.end = end
        self.error_score = error_score
        self.coverage = end - start

    def __repr__(self) -> str:
        return (
            f"Range({self.id}, {self.start}, {self.end}, "
            f"{self.error_score}, coverage = {self.coverage})"
        )


def select_ranges(ranges: list) -> list:
    """
    Select a non-overlapping set of ranges.

    This is used to choose PDB chains that cover the protein sequence while
    preferring better error scores and longer coverage.

    Parameters
    ----------
    ranges : list
        List of Range objects.

    Returns
    -------
    list
        Selected ranges.

    Examples
    --------
    >>> rs = [Range("A", 1, 10, 3.0), Range("B", 5, 12, 2.0)]
    >>> out = select_ranges(rs)
    >>> isinstance(out, list)
    True
    """
    sorted_ranges = sorted(ranges, key=lambda r: (r.start, r.end))
    selected_ranges: List[Range] = []
    last_end = float("-inf")

    for current_range in sorted_ranges:
        if current_range.start > last_end:
            selected_ranges.append(current_range)
            last_end = current_range.end
        else:
            if (
                current_range.error_score < selected_ranges[-1].error_score
                or abs(current_range.error_score - selected_ranges[-1].error_score) < 0.5
            ):
                if current_range.coverage > selected_ranges[-1].coverage:
                    selected_ranges.append(current_range)
                last_end = max(last_end, current_range.end)
    return selected_ranges


def get_optimal_pdbs_df(
    df: pd.DataFrame, seq_length: int, min_res_val: float = 2.5
) -> pd.DataFrame:
    """
    Choose a subset of PDB structures to maximize coverage with good resolution.

    Parameters
    ----------
    df : pandas.DataFrame
        Input PDB metadata table.
    seq_length : int
        Protein sequence length. Used to derive coverage comparisons.
    min_res_val : float, default 2.5
        Resolution cutoff used for initial filtering.

    Returns
    -------
    pandas.DataFrame
        Subset of the input dataframe containing selected PDB rows.

    Examples
    --------
    >>> out = get_optimal_pdbs_df(pd.DataFrame(), seq_length=100)
    """
    df = df.copy()
    df["range_obj"] = df.apply(
        lambda x: Range(x.id, x.chain_start, x.chain_end, x.resolution), axis=1
    )
    min_res_for_max_cov = min(df[df["coverage"] == max(df["coverage"])]["resolution"].values)

    df_1 = df[df["resolution"] <= min_res_val]
    df_2 = df[df["resolution"] <= min_res_for_max_cov]

    data1 = df_1["range_obj"].to_list()
    data2 = df_2["range_obj"].to_list()

    selected_ranges1 = select_ranges(data1)
    selected_ranges2 = select_ranges(data2)

    selected_ranges = list(set(selected_ranges1 + selected_ranges2))
    ids = [i.id for i in selected_ranges]
    selected_pdbs = df.loc[ids]
    return selected_pdbs


def download_pdbs(gene_name: str, gene_dir: str, pdb_id_list: list) -> list:
    """
    Download PDB files from PDBe into a gene cache directory.

    Parameters
    ----------
    gene_name : str
        Gene symbol used in the filename prefix.
    gene_dir : str
        Output directory where PDB files will be written.
    pdb_id_list : list
        List of PDB identifiers to download.

    Returns
    -------
    list
        PDB ids that failed to download.

    Examples
    --------
    >>> failed = download_pdbs("TP53", "/tmp/tp53", ["1ABC"])
    """
    failed_pdbs: List[str] = []
    os.makedirs(gene_dir, exist_ok=True)
    for pdb_id in pdb_id_list:
        target = os.path.join(gene_dir, f"g_{gene_name}_p_{pdb_id}.pdb")
        if os.path.isfile(target):
            continue
        url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pdb_id.lower()}.ent"
        try:
            resp = requests.get(url, timeout=60)
            if resp.status_code != 200 or not resp.content:
                raise RuntimeError(f"bad response status {resp.status_code} for {pdb_id}")
            with open(target, "wb") as f:
                f.write(resp.content)
        except Exception as e:  # noqa: BLE001
            if os.path.isfile(target):
                os.remove(target)
            failed_pdbs.append(pdb_id)
            log_progress("pdb_clean", 0, f"download failed for {pdb_id}: {e}")
    return failed_pdbs


def get_chain_ids(pdb_file: str) -> list:
    """
    Extract chain identifiers from a PDB file.

    Parameters
    ----------
    pdb_file : str
        Path to a PDB file.

    Returns
    -------
    list
        Chain ids found in the structure.

    Examples
    --------
    >>> ids = get_chain_ids("protein.pdb")
    """
    try:
        from Bio import PDB  # type: ignore
    except ImportError:
        raise ImportError("Biopython is required for pdb_clean but not installed")
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("PDB_structure", pdb_file)
    chain_ids = set()
    for model in structure:
        for chain in model:
            chain_ids.add(chain.id)
    return list(chain_ids)


def _clean_gene_pdbs(
    gene_name: str,
    gene_dir: str,
    id_chain_map: Dict[str, Tuple[List[str], List[str]]],
) -> Tuple[List[Optional[str]], List[str]]:
    """
    Clean downloaded PDBs for a gene.

    Parameters
    ----------
    gene_name : str
        Gene symbol used in the filename prefix.
    gene_dir : str
        Directory containing downloaded PDB files.
    id_chain_map : dict[str, tuple[list[str], list[str]]]
        Mapping of PDB id to (kept_chains, removed_chains).

    Returns
    -------
    tuple[list[str or None], list[str]]
        Cleaned paths (or None for failures) and the list of PDB ids that could
        not be cleaned.

    Examples
    --------
    >>> paths, bad = _clean_gene_pdbs("TP53", "/tmp/tp53", {})
    """
    cleaned_pdbs_path: List[Optional[str]] = []
    uncleanable_pdbs: List[str] = []
    for pdb_id, chains in id_chain_map.items():
        input_path = os.path.join(gene_dir, f"g_{gene_name}_p_{pdb_id}.pdb")
        output_path = os.path.join(gene_dir, f"cleaned_g_{gene_name}_p_{pdb_id}.pdb")
        out = clean_pdb_file(
            input_pdb_path=input_path,
            output_pdb_path=output_path,
            remove_chains=chains[-1],
            remove_heterogens=True,
            remove_water=True,
            add_hydrogens=True,
            ph=7.0,
        )
        cleaned_pdbs_path.append(out)
        if out is None:
            uncleanable_pdbs.append(pdb_id)
    return cleaned_pdbs_path, uncleanable_pdbs


def _run_pdb_clean_for_gene(
    gene_name: str,
    entry_id: str,
    gene_dir: str,
    min_res_val: float = 2.5,
) -> pd.DataFrame:
    """
    Download, select, and clean PDB structures for one gene.

    Parameters
    ----------
    gene_name : str
        Gene symbol.
    entry_id : str
        UniProt entry id.
    gene_dir : str
        Output directory for all per-gene cached files.
    min_res_val : float, default 2.5
        Resolution cutoff used during selection.

    Returns
    -------
    pandas.DataFrame
        Selected PDB metadata with a local path column.

    Examples
    --------
    >>> df = _run_pdb_clean_for_gene("GBA3", "Q9H3H0", "/tmp/gba3")
    """
    max_coverage, _prot_seq = get_protein_details(entry_id)
    pdbs_df = get_pdbs_df(entry_id)

    failed_pdbs: Optional[List[str]] = None
    selected_pdbs: Optional[pd.DataFrame] = None
    uncleanable_pdbs: List[str] = []
    cleaned_pdbs_path: List[Optional[str]] = []

    while failed_pdbs != []:
        if failed_pdbs is not None:
            pdbs_df = pdbs_df.drop(index=failed_pdbs, errors="ignore")
        if len(pdbs_df) == 0:
            raise OverflowError("no pdbs left to check, likely over sized pdbs")
        selected_pdbs = get_optimal_pdbs_df(pdbs_df, max_coverage, min_res_val=min_res_val)
        id_chain_map_raw: Dict[str, List[str]] = selected_pdbs["chain_type"].to_dict()
        failed_pdbs = download_pdbs(gene_name, gene_dir, list(id_chain_map_raw.keys()))
        if failed_pdbs:
            continue

        id_chain_map: Dict[str, Tuple[List[str], List[str]]] = {}
        for pdb_id, chains in id_chain_map_raw.items():
            all_chains = get_chain_ids(os.path.join(gene_dir, f"g_{gene_name}_p_{pdb_id}.pdb"))
            chains_to_remove = list(set(all_chains).difference(set(chains)))
            id_chain_map[pdb_id] = (chains, chains_to_remove)

        cleaned_pdbs_path, uncleanable_pdbs = _clean_gene_pdbs(gene_name, gene_dir, id_chain_map)
        failed_pdbs = list(uncleanable_pdbs)

    if selected_pdbs is None:
        raise RuntimeError("Failed to select optimal PDBs for gene {0}".format(gene_name))

    if not cleaned_pdbs_path:
        raise Exception("No clean pdbs")

    selected_pdbs = selected_pdbs.copy()
    selected_pdbs["path"] = cleaned_pdbs_path
    selected_pdbs.to_csv(os.path.join(gene_dir, f"{gene_name}_pdbs.csv"))
    return selected_pdbs


def _upload_gene_artifacts(
    gene_name: str,
    scan_id: str,
    output: str,
    gene_dir: str,
    selected_pdbs: pd.DataFrame,
) -> Dict[str, object]:
    """
    Upload per-gene artifacts to the datastore.

    Parameters
    ----------
    gene_name : str
        Gene symbol.
    scan_id : str
        Scan identifier used to group outputs.
    output : str
        Output name prefix used for uploaded artifacts.
    gene_dir : str
        Local cache directory for this gene.
    selected_pdbs : pandas.DataFrame
        Selected PDB rows with a local path column.

    Returns
    -------
    dict
        Metadata describing uploaded addresses and local paths.

    Examples
    --------
    >>> meta = _upload_gene_artifacts("TP53", "scan1", "out", "/tmp/tp53", pd.DataFrame())
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")

    pdb_addresses: Dict[str, str] = {}
    for pdb_id in list(selected_pdbs.index):
        cleaned_path = selected_pdbs.loc[pdb_id]["path"]
        if cleaned_path is None or not os.path.exists(cleaned_path):
            continue
        with open(cleaned_path, "r") as f:
            cleaned_content = f.read()
        card = DataCard(address="", file_type="pdb", data_type="text/plain")
        pdb_key = f"{output}_cleaned_{gene_name}_{pdb_id}.pdb"
        pdb_addresses[pdb_id] = datastore.upload_data_from_memory(cleaned_content, pdb_key, card)

    csv_card = DataCard(address="", file_type="csv", data_type="pandas.DataFrame")
    csv_key = f"{output}_{gene_name}_pdbs.csv"
    csv_address = datastore.upload_data_from_memory(selected_pdbs, csv_key, csv_card)

    metadata = {
        "gene_name": gene_name,
        "scan_id": scan_id,
        "local_gene_dir": gene_dir,
        "pdbs_csv_address": csv_address,
        "num_pdbs": int(len(selected_pdbs)),
        "selected_pdb_ids": [str(i) for i in list(selected_pdbs.index)],
        "pdb_addresses": pdb_addresses,
        "records": selected_pdbs.reset_index(drop=True).to_dict(orient="records"),
    }
    return metadata


def pdb_clean(
    gene_name: str,
    entry_id: str,
    scan_id: str,
    output: str,
    min_res_val: float = 2.5,
) -> str:
    """Fetch and clean PDB structures for a gene.

    Parameters
    ----------
    gene_name: str
        Gene symbol, for example "GBA3".
    entry_id: str
        UniProt entry id, for example "Q9H3H0".
    scan_id: str
        Logical id that groups all files produced by a single proteome
        scan run. Becomes the top-level folder under the cache root.
    output: str
        Output name prefix for uploaded datastore artifacts.
    min_res_val: float
        Minimum resolution cutoff used during selection (default 2.5).

    Returns
    -------
    str
        DeepChem address of the metadata JSON, which references the
        per-gene CSV and per-PDB cleaned file addresses.

    Examples
    --------
    >>> addr = pdb_clean("GBA3", "Q9H3H0", "scan1", "my_run")
    """
    datastore = config.get_datastore()
    if datastore is None:
        raise ValueError("Datastore not set")
    if not gene_name:
        raise ValueError("gene_name is required")
    if not entry_id:
        raise ValueError("entry_id is required")
    if not scan_id:
        raise ValueError("scan_id is required")
    if not output:
        raise ValueError("output is required")

    gene_dir = str(ps_cache.get_gene_dir(scan_id, gene_name))
    csv_path = os.path.join(gene_dir, f"{gene_name}_pdbs.csv")

    log_progress("pdb_clean", 5, f"cache dir for {gene_name}: {gene_dir}")

    if os.path.exists(csv_path):
        log_progress("pdb_clean", 20, f"resuming from existing {csv_path}")
        selected_pdbs = pd.read_csv(csv_path, index_col="id")
    else:
        log_progress("pdb_clean", 20, f"fetching PDBs for {gene_name} ({entry_id})")
        selected_pdbs = _run_pdb_clean_for_gene(
            gene_name, entry_id, gene_dir, min_res_val=min_res_val
        )

    log_progress("pdb_clean", 70, "uploading gene artifacts to datastore")
    metadata = _upload_gene_artifacts(
        gene_name=gene_name,
        scan_id=scan_id,
        output=output,
        gene_dir=gene_dir,
        selected_pdbs=selected_pdbs,
    )

    metadata_json = json.dumps(metadata, indent=2, default=str)
    json_card = DataCard(address="", file_type="json", data_type="json")
    metadata_key = f"{output}_{gene_name}_metadata.json"
    metadata_address = datastore.upload_data_from_memory(metadata_json, metadata_key, json_card)

    log_progress("pdb_clean", 100, "pdb_clean complete")
    return metadata_address


__all__ = [
    "pdb_clean",
    "get_protein_details",
    "get_pdbs_df",
    "get_optimal_pdbs_df",
    "select_ranges",
    "Range",
    "download_pdbs",
    "get_chain_ids",
    "clean_pdb_file",
]
