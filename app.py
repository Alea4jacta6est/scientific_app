import py3Dmol
import streamlit as st
from stmol import showmol
import numpy as np
from protein_utils import (
    count_distance,
    get_amino_acids,
    get_mhc_chains,
    get_prosite_positions,
    get_pseudosequence,
)

cartoon_radius = 0.2
stick_radius = 0.2
styles = ["cartoon", "line", "cross", "stick", "sphere"]


def main():
    st.set_page_config(page_title="PDBViz", layout="wide", initial_sidebar_state="auto")
    st.title("Protein data bank visualization")
    pdb_id = st.sidebar.text_area("Enter PDB id", "1A4J").strip()
    pattern_id = st.sidebar.text_area(f"Enter PROSITE pattern id", "PS00290").strip()
    threshold = st.sidebar.number_input("Distance threshold", 10)
    bcolor = st.sidebar.color_picker("Pick a background color", "#FFFFFF")
    style = st.sidebar.selectbox("Select a style", styles)
    chosen_color = st.sidebar.text_input("Select a color for PROSITE pattern", "red")
    try:
        sequence = get_amino_acids(pdb_id)[0]
        chosen_chain = st.sidebar.selectbox(
            "Select a chain to highlight", [c.id for c in get_mhc_chains(pdb_id)]
        )
        view = py3Dmol.view(query="pdb:" + pdb_id)
        view.setStyle({style: {"color": "white"}})
        view.setBackgroundColor(bcolor)
    except ValueError:
        st.error("You provided wrong pdb id")

    try:
        positions = get_prosite_positions(pattern_id, sequence)
        for hl_resi in positions:
            view.addStyle(
                {"chain": chosen_chain, "resi": hl_resi, "elem": "C"},
                {style: {"color": chosen_color, "radius": stick_radius}},
            )

            view.addStyle(
                {"chain": chosen_chain, "resi": hl_resi},
                {style: {"radius": stick_radius}},
            )
        st.subheader(f"3D structure of a PDB file")
        showmol(view, width=1200)
        pseudoseq = get_pseudosequence(sequence, positions)
        st.subheader(f"Pseudosequence")
        st.text(pseudoseq)
        if len(positions) > 1:
            st.text(f"Start: {positions[0]} Stop: {positions[-1]}")
        st.subheader("A matrix with distances")
        chains_df = count_distance(get_mhc_chains(pdb_id))
        subset_df = (
            chains_df.where(chains_df >= threshold, np.nan)
            .dropna(how="all", axis=0)
            .dropna(how="all", axis=1)
        )
        st.dataframe(subset_df)
    except:
        st.error("You provided wrong prosite id")


if __name__ == "__main__":
    main()
