from .calculators import *

""" 
TO-DO:
    - [x] Add distance
    - [x] Add angle
    - [x] Add dihedral
    - [x] Add planar_angle
    - [x] Add contacts
    - [x] Add RMSD
    - [x] Add distWATbridge
    - [x] Add pKa
"""


def run_distance(Measure):
    """
    DESCRIPTION:

    """

    Measure.result.append(
        calc_distance(Measure.sel[0], Measure.sel[1], Measure.options["type"])
    )


def run_angle(Measure):
    """
    DESCRIPTION:

    """

    Measure.result.append(
        calc_angle(
            Measure.sel[0],
            Measure.sel[1],
            Measure.sel[2],
            Measure.options["units"],
            Measure.options["domain"],
        )
    )


def run_dihedral(Measure):
    """
    DESCRIPTION:

    """

    Measure.result.append(
        calc_dihedral(
            Measure.sel[0],
            Measure.sel[1],
            Measure.sel[2],
            Measure.sel[3],
            Measure.options["units"],
            Measure.options["domain"],
        )
    )


def run_planar_angle(Measure):
    """
    DESCRIPTION:

    """

    Measure.result.append(
        calc_planar_angle(
            Measure.sel[0],
            Measure.sel[1],
            Measure.options["units"],
            Measure.options["domain"],
        )
    )


def run_pka(Measure):
    """
    DESCRIPTION:

    """

    Measure.result.append(
        calc_pka(
            Measure.sel[0],
            Measure.options["pka_ref"],
            Measure.options["pdb_folder"],
            Measure.options["keep_pdb"],
            Measure.options["keep_pka"],
        )
    )


def run_contacts(Measure):
    """
    DESCRIPTION:

    """

    if Measure.options["mode"] == "selection":
        Measure.result.append(
            calc_contacts_selection(
                Measure.sel[0],
                Measure.sel[1],
                Measure.options["interactions"],
                Measure.options["measure_dists"],
                Measure.options["out_format"],
            )
        )

    elif Measure.options["mode"] == "protein":
        Measure.result.append(
            calc_contacts_protein(
                Measure.sel[0],
                Measure.sel[1],
                Measure.options["interactions"],
                Measure.options["measure_dists"],
                Measure.options["out_format"],
            )
        )


def run_RMSD(Measure):
    """
    DESCRIPTION:

    """

    Measure.result.append(
        calc_RMSD(
            Measure.sel[0], Measure.options["ref"], Measure.options["superposition"]
        )
    )


def run_distWATbridge(Measure):
    """
    DESCRIPTION:

    """

    Measure.result.append(
        calc_distWATbridge(
            Measure.sel[0],
            Measure.sel[1],
            Measure.sel[2],
            Measure.sel[3],
            Measure.sel[4],
            Measure.sel[5],
        )
    )
