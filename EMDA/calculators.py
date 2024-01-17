from MDAnalysis.core.groups import AtomGroup
import MDAnalysis.lib.distances as mdadist
import MDAnalysis.analysis.rms as rms

from ..exceptions import (
    NotExistingInteraction,
)

from ..snippets import check_folder

from subprocess import run as run_command
from os import remove, chdir


"""
DESCRIPTION
    This Python file contains the calculators that are called from measurements.py.
    Each function needs to input the selections and the options.
"""

def build_plane(
    positions,  # nx3, being n the number of atoms that build the plane
):
    """
    DESCRIPTION
        Function for obtaining the equation of the plane from a given list containing the positions of three atoms.
        It is used in the 'planar_ange' function.
    """
    from numpy import cross, dot

    if len(positions) == 3:
        v1 = positions[2] - positions[0]
        v2 = positions[1] - positions[0]

        a, b, c = cross(v1, v2)
        d = dot([a, b, c], positions[2])

        #        return array(a, b, c, d)
        return [a, b, c, d]


def planar_angle(
    plane_A,  # nx3 np array, being n the number of atoms for building the plane (min 3)
    plane_B,  # nx3 np array, being n the number of atoms for building the plane (min 3)
    units,  # [ deg | rad ]
    domain,  # [ 360 | 180 ] 360: (0, 360); 180: (-180, 180)
):
    from numpy import arccos, sqrt

    plane_A = build_plane(plane_A)
    plane_B = build_plane(plane_B)

    ang = arccos(
        (plane_A[0] * plane_B[0] + plane_A[1] * plane_B[1] + plane_A[2] * plane_B[2])
        / (
            sqrt(plane_A[0] ** 2 + plane_A[1] ** 2 + plane_A[2] ** 2)
            * sqrt(plane_B[0] ** 2 + plane_B[1] ** 2 + plane_B[2] ** 2)
        )
    )

    if units in ("rad", "radian", "radians"):

        if domain in (180, "180", "pi"):
            pass

        elif domain in (360, "360", "2pi"):
            from math import pi

            ang = (ang + pi) % pi

    elif units in ("deg", "degree", "degrees"):
        from numpy import rad2deg

        if domain in (180, "180", "pi"):
            ang = rad2deg(ang)

        elif domain in (360, "360", "2pi"):
            ang = (rad2deg(ang) + 360) % 360

    return ang


def distance(
    sel1, # selection containing one or more atoms
    sel2, # selection containing one or more atoms
    type # [min | max | com | cog ]
):
    """
    DESCRIPTION
        Function that takes two selections and gives:
            - minimum distance
            - maximum distance
            - distance between Centers of Mass (com)
            - distance between Centers of Geometry (cog)
    """
    from numpy import min as npmin
    from numpy import max as npmax

    if type == 'min':
        d = npmin(mdadist.distance_array(sel1, sel2, backend='OpenMP'))

    elif type == 'max':
        d = npmax(mdadist.distance_array(sel1, sel2, backend='OpenMP'))

    elif type == 'com':
        d = mdadist.distance_array(sel1.center_of_mass(), sel2.center_of_mass(), backend='OpenMP')

    elif type == 'cog':
        d = mdadist.distance_array(sel1.center_of_geometry(), sel2.center_of_geometry(), backend='OpenMP')

    return d


def dihedral(
    sel1,
    sel2,
    sel3,
    sel4,
    units, # [ rad | deg ]
    domain # [ 180 | 360 ]
):
    """
    DESCRIPTION
        Function that calculates the dihedral angle between four atoms.
        Dihedral angle is calculated by measuring the planar angle of the planes builds by the sel1, sel2, sel3 and sel2, sel3, sel4.

    OPTIONS (as arguments)
        - units: radians (rad) or degrees (deg)
        - domain: -180 to 180 º or 0 to 360º
    """

    d = mdadist.calc_dihedrals(sel1, sel2, sel3, sel4, backend="OpenMP")

    if units in ("rad", "radian", "radians", "pi"):
        from math import pi

        if domain in (360, "360", "2pi"):
            d = ( d + pi ) % pi
            return d
        elif domain in (180, "180", "pi"):
            return d

    elif units in ("deg", "degree", "degrees"):
        from numpy import rad2deg

        d = rad2deg(d)
        if domain in (360, "360", "2pi"):
            d = ( d + 360 ) % 360
            return d
        elif domain in (180, "180", "pi"):
            return d


def angle(
    sel1,
    sel2,
    sel3,
    units, # [ rad | deg ]
    domain # [ 180 | 360 ]
):
    """
    DESCRIPTION
        Function that calculates the angle between three atoms.
        Order of selection is important (sel2 is the vertex)

    OPTIONS (as arguments)
        - units: radians (rad) or degrees (deg)
        - domain: -180 to 180 º or 0 to 360º
    """

    a = mdadist.calc_angles(sel1, sel2, sel3, backend="OpenMP")

    if units in ("rad", "radian", "radians", "pi"):
        from math import pi

        if domain in (360, "360", "2pi"):
            a = ( a + pi ) % pi
            return a
        elif domain in (180, "180", "pi"):
            return a

    elif units in ("deg", "degree", "degrees"):
        from numpy import rad2deg

        a = rad2deg(a)
        if domain in (360, "360", "2pi"):
            a = ( a + 360 ) % 360
            return a
        elif domain in (180, "180", "pi"):
            return a
        

def pka(
        sel_protein,            # whole protein AtomGroup for saving PDB
        #sel_pka=None,           # optional selection for checking only the pKa of the selected resids
        pka_ref='neutral',  
        pdb_folder='.propka',  # folder to store files
        frame='current',
        keep_pdb=False,
        keep_pka=False,
):
    

    # save pdb
    chdir(pdb_folder)
    sel_protein.write(f"{frame}.pdb")

    # predict pKa with PROpKa
    run_command(
        [
            "propka3",
            "-r", pka_ref,
            f"{frame}.pdb",
        ]
    )

    # read .pka file
    f = open(f"{frame}.pka").readlines()

    # parse .pka file
    save_pka = False
    for l in f:
        if l.startswith('SUMMARY OF THIS PREDICTION'):
            save_pka = True
            pkas = {}
            continue

        if save_pka:
            if l.startswith('--------------------------------------------------------------------------------------------------------'):
                save_pka = False
                continue
            elif l.startswith('       Group'):
                pass
            else :
                l = l.split()
                if float(l[3]) != 99.99:
                    pkas[l[0] + str(int(l[1]))] = float(l[3])

    # remove PDBs if indicated
    if not keep_pdb:
        remove(f"{frame}.pdb")

    # remove .pKas if indicated
    if not keep_pka:
        remove(f"{frame}.pka")
    
    chdir ('..')

    return pkas


def contacts_selection(
        sel,
        sel_env,
        interactions,
        measure_distances=False,
        out_format='new',
):

    residues = []
    for name, id in zip(sel_env.resnames, sel_env.resids):
        if {"name" : name, "id" : id} not in residues:
            residues.append({"name" : name, "id" : id})

    contacts = {}

    for r in range(len(residues)):
        if residues[r]["id"] not in sel.residues.resindices: #residue[1] is the residue's number
            if residues[r]["name"] in interactions:
                if out_format == 'old':
                    contacts[residues[r]['id']] = residues[r]['name']
                
                elif out_format == 'new':
                    if measure_distances:
                        contacts[residues[r]['name'] + str(residues[r]['id'])] = distance(sel, sel_env.residues[r].atoms, type='min')
                    elif not measure_distances:
                        contacts[residues[r]['name'] + str(residues[r]['id'])] = None

    return contacts


def contacts_protein(
        sel, 
        sel_env,
        interactions,
        measure_distances= False,
        out_format= 'new'
):
    
    contacts = {}

    protein_residues = ["ARG", "HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "LYS", "ASP", "ASH", "GLU", "GLH", "SER", "THR", 
                        "ASN", "GLN", "CYS", "SEC", "GLY", "PRO", "ALA", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "VAL"]
    

    for residue in sel.atoms.select_atoms('protein').residues:
        #if residue.resname not in protein_residues:
            # if the residue does not belong to the protein, skip it
        #    continue

        residue_env = sel.atoms.select_atoms('around %s group select' % str(sel_env), select=residue.atoms)

        contacts[residue.resname + str(residue.resid)] = contacts_selection(
            sel=residue.atoms,
            sel_env=residue_env,
            interactions=interactions,
            measure_distances=measure_distances,
            out_format=out_format
        )

    return contacts




