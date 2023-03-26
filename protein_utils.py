import re
from typing import List

import pandas as pd
import requests
from Bio.ExPASy import Prosite, get_prosite_raw
from Bio.PDB import PDBList, PDBParser, is_aa
from scipy.spatial.distance import euclidean


def get_amino_acids(code: str):
    link = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{code}"
    try:
        data = requests.get(link).json()[code.lower()]
        return [data[i]["sequence"] for i in range(len(data)) if "sequence" in data[i]]
    except Exception:
        raise ValueError("Wrong pdb id, double check your input")


def get_mhc_chains(pdb_id):
    parser = PDBParser()
    pdb_file = PDBList().retrieve_pdb_file(pdb_id, pdir="./files", file_format="pdb")
    structure = parser.get_structure(pdb_id, pdb_file)
    return structure.get_chains()


def get_prosite_positions(pattern_id: str, protein_sequence: str):
    handle = get_prosite_raw(pattern_id)
    record = Prosite.read(handle)
    prosite_pattern = record.pattern
    prosite_pattern = convert_prosite_to_regex(prosite_pattern)
    matches = re.finditer(prosite_pattern, protein_sequence)
    highlighted_positions = []
    for match in matches:
        for i in range(match.start(), match.end()):
            highlighted_positions.append(i)
    return highlighted_positions


def convert_prosite_to_regex(prosite_pattern: str):
    # Replace PROSITE residue codes with regular expression character classes
    regex_pattern = prosite_pattern.replace("X", "[A-Z]")
    regex_pattern = regex_pattern.replace("x", "[a-z]")
    regex_pattern = regex_pattern.replace("Z", "[EQ]")
    regex_pattern = regex_pattern.replace("B", "[DN]")
    regex_pattern = regex_pattern.replace("J", "[IL]")
    regex_pattern = regex_pattern.replace("O", "[KQRN]")
    regex_pattern = regex_pattern.replace("U", "[CYS]")
    # Replace PROSITE special characters with regular expression syntax
    regex_pattern = regex_pattern.replace("{", "[^")
    regex_pattern = regex_pattern.replace("}", "]")
    regex_pattern = regex_pattern.replace("-", "")
    return regex_pattern


def get_pseudosequence(seq: str, mhc_residues: List[int]) -> str:
    """Transforms the original sequence to a pseudosequence

    Args:
        seq (str): original sequence
        mhc_residues (list): positions nums

    Returns:
        pseudoseq (str): _description_
    """
    pseudoseq = ""
    for i, residue in enumerate(seq):
        if i + 1 not in mhc_residues:
            pseudoseq += "-"
        else:
            pseudoseq += residue
    return pseudoseq


def count_distance(mhc_chains):
    chain_L, chain_H = None, None
    for chain in mhc_chains:
        if chain.id == "L":
            chain_L = chain
        elif chain.id == "H":
            chain_H = chain
    if chain_L is None or chain_H is None:
        raise ValueError("Could not find L and H chains in MHC chains")
    distances = pd.DataFrame(
        index=[
            residue.id[1]
            for residue in chain_L.get_residues()
            if is_aa(residue.get_resname(), standard=True)
        ],
        columns=[
            residue.id[1]
            for residue in chain_H.get_residues()
            if is_aa(residue.get_resname(), standard=True)
        ],
    )
    for residue_L in chain_L.get_residues():
        if not is_aa(residue_L.get_resname(), standard=True):
            continue
        for residue_H in chain_H.get_residues():
            if not is_aa(residue_H.get_resname(), standard=True):
                continue
            distance = euclidean(
                residue_L["CA"].get_coord(), residue_H["CA"].get_coord()
            )
            distances.at[residue_L.id[1], residue_H.id[1]] = distance
    distances = distances.astype(float)
    return distances


if __name__ == "__main__":
    pdb_id = "1A4J"
    sequence = get_amino_acids(pdb_id)[0]
    mhc_chains = get_mhc_chains(pdb_id)
    print(count_distance(mhc_chains))
