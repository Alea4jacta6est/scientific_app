import streamlit as st

from viz_protein import show_st_3dmol


def main():
    st.set_page_config(page_title="PDBViz", layout="wide", initial_sidebar_state="auto")
    st.title("Protein data bank visualization")
    pdb_id = st.sidebar.text_area("Enter PDB id", "4HHB").strip()
    st.subheader(f"3D structure of a PDB file")
    show_st_3dmol(pdb_id)

    prosite = st.sidebar.text_area(f"Enter PROSITE pattern id", "").strip()
    dist = st.sidebar.number_input("Distance threshold", 0.1)


if __name__ == "__main__":
    main()
