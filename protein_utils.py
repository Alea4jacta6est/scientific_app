import requests
import re
import numpy as np
from Bio.ExPASy import Prosite, get_prosite_raw
from Bio.PDB import PDBParser, PDBList
from typing import List


def get_amino_acids(code: str):
    try:
        data = requests.get(
            f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{code}"
        ).json()[code.lower()]
        return [data[i]["sequence"] for i in range(len(data)) if "sequence" in data[i]]
    except Exception:
        raise ValueError("Wrong pdb id, double check your input")


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
        if i + 1 in mhc_residues:
            pseudoseq += "-"
        else:
            pseudoseq += residue
    return pseudoseq


# def get_filtered_residues(pdb_id, dist_threshold):
#     d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
#     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
#     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
#     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
#     parser = PDBParser()
#     pdbl = PDBList()
#     pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir='./files', file_format='pdb')
#     peptide_struct = parser.get_structure(pdb_id, pdb_file)
#     for model in peptide_struct:
#          for chain in model:
#             seq = []
#             for residue in chain:
#                 if residue.resname in d3to1:
#                     seq.append(d3to1[residue.resname])
#     print("".join(seq))
#     # Extract peptide and MHC chain objects
#     peptide = peptide_struct[0]['A']
#     print([i for i in peptide.get_residues()])
# print(dir(peptide))
# # # mhc_chain = mhc_struct[0]['A']
# # Get coordinates of all atoms in the peptide and MHC chain
# peptide_atoms = [atom.get_coord() for atom in peptide.get_atoms() if atom.get_name() == 'CA']
# mhc_atoms = [atom.get_coord() for atom in mhc_chain.get_atoms() if atom.get_name() == 'CA']
# # Calculate pairwise distances between all atoms in the peptide and MHC chain
# dist_matrix = np.linalg.norm(np.array(peptide_atoms)[:, np.newaxis, :] - np.array(mhc_atoms)[np.newaxis, :, :], axis=-1)

# # Get minimum distances for each amino acid in the peptide
# min_dists = np.min(dist_matrix, axis=1)
# filtered_residues = [residue for residue, dist in zip(peptide, min_dists) if dist <= dist_threshold]
# return filtered_residues

# print(get_filtered_residues(pdb_id="1A4J", dist_threshold=0.3))
