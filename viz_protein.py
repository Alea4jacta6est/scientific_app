import py3Dmol
from stmol import showmol


def show_st_3dmol(
    pdb_code,
    style_lst=None,
    label_lst=None,
    reslabel_lst=None,
    zoom_dict=None,
    surface_lst=None,
    cartoon_style="trace",
    cartoon_radius=0.2,
    cartoon_color="lightgray",
    zoom=1,
    spin_on=False,
    width=900,
    height=600,
):

    """
    Show 3D view of protein structure from the
    Protein Data Bank (PDB)

    Parameters
    ----------
    pdb_code: str
        Four-letter code of protein structure in the PDB
        (e.g., 5P21)
    style_lst: list of lists of dicts
        A nested list with each sublist containing a
        selection dictionary at index 0 and coloring
        dictionary at index 1
    label_lst: list of lists of dicts
        A nested list with each sublist containing a
        label string at index 0, coloring dictionary
        at index 1, and selection dictionary at
        index 2
    reslabel_lst: list of lists of dicts
        A nested list with each sublist containing a
        selection dictionary at index 0 and coloring
        dictionary at index 1
    zoom_dict: dict
    surface_lst: list of lists of dicts
        A nested list with each sublist containing a
        coloring dictionary at index 0 and selection
        dictionary at index 1
    cartoon_style: str
        Style of protein structure backbone cartoon
        rendering, which can be "trace", "oval", "rectangle",
        "parabola", or "edged"
    cartoon_radius: float
        Radius of backbone cartoon rendering
    cartoon_color: str
        Color of backbone cartoon rendering
    zoom: float
        Level of zoom into protein structure
        in unit of Angstroms
    spin_on: bool
        Boolean specifying whether the visualized
        protein structure should be continually
        spinning (True) or not (False)
    width: int
        Width of molecular viewer
    height: int
        Height of molecular viewer
    """
    view = py3Dmol.view(query=f"pdb:{pdb_code.lower()}", width=width, height=height)

    view.setStyle(
        {
            "cartoon": {
                "style": cartoon_style,
                "color": cartoon_color,
                "thickness": cartoon_radius,
            }
        }
    )

    if surface_lst is not None:
        for surface in surface_lst:
            view.addSurface(py3Dmol.VDW, surface[0], surface[1])

    if style_lst is not None:
        for style in style_lst:
            view.addStyle(
                style[0],
                style[1],
            )

    if label_lst is not None:
        for label in label_lst:
            view.addLabel(label[0], label[1], label[2])

    if reslabel_lst is not None:
        for reslabel in reslabel_lst:
            view.addResLabels(reslabel[0], reslabel[1])

    if zoom_dict is None:
        view.zoomTo()
    else:
        view.zoomTo(zoom_dict)

    view.spin(spin_on)

    view.zoom(zoom)
    showmol(view, height=height, width=width)
