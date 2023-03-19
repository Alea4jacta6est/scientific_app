import streamlit as st
import py3Dmol
from stmol import showmol
from protein_utils import get_amino_acids, get_prosite_positions, get_pseudosequence


cartoon_radius = 0.2
stick_radius = 0.2


def main():
    st.set_page_config(page_title="PDBViz", layout="wide", initial_sidebar_state="auto")
    st.title("Protein data bank visualization")
    pdb_id = st.sidebar.text_area("Enter PDB id", "1A4J").strip()
    pattern_id = st.sidebar.text_area(f"Enter PROSITE pattern id", "PS00290").strip()
    dist = st.sidebar.number_input("Distance threshold", 0.3)
    bcolor = st.sidebar.color_picker("Pick A Color", "#FFFFFF")
    style = st.sidebar.selectbox(
        "style", ["cartoon", "line", "cross", "stick", "sphere"]
    )
    hl_color = st.sidebar.text_input(label="Highlight Color", value="red")
    hl_chain = st.sidebar.text_input(label="Highlight Chain", value="A")
    try:
        sequence = get_amino_acids(pdb_id)[0]
        view = py3Dmol.view(query="pdb:" + pdb_id)
        view.setStyle({style: {"color": "white"}})
        view.setBackgroundColor(bcolor)
    except ValueError:
        st.error("You provided wrong pdb id")

    try:
        positions = get_prosite_positions(pattern_id, sequence)
        for hl_resi in positions:
            view.addStyle(
                {"chain": hl_chain, "resi": hl_resi, "elem": "C"},
                {style: {"color": hl_color, "radius": stick_radius}},
            )

            view.addStyle(
                {"chain": hl_chain, "resi": hl_resi}, {style: {"radius": stick_radius}}
            )
        st.subheader(f"3D structure of a PDB file")
        showmol(view, width=1200)
        pseudoseq = get_pseudosequence(sequence, positions)
        st.subheader(f"Pseudosequence")
        st.text(pseudoseq)
        if len(positions) > 1:
            st.text(f"Start: {positions[0]} Stop: {positions[-1]}")
    except:
        st.error("You provided wrong prosite id")


if __name__ == "__main__":
    main()
