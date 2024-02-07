# from dataclasses import dataclass

from MDAnalysis.core.groups import AtomGroup

# load custom exceptions
from .exceptions import NotAvailableOptionError
from .exceptions import NotSingleAtomSelectionError
from .exceptions import NotThreeAtomsSelectionError
from .exceptions import NotExistingInteractionError
from .exceptions import NotExistingSelectionError

from .selection import convert_selection

# @dataclass
# class Measure:
#    name    : str
#    type    : str
#    sel     : list
#    options : dict
#    result  : list


"""
HOW TO BUILD AN ADDER:

    - name and sel(s) have to be requested, other options can be added too
    - Check that all the input data is in accordance with the created calculator
    - finish the adder by adding the configuration as a Measure data class in the measurers list of the EMDA class like so:
        
        self.measures[name] = self.Measure(
            name    = name,
            type    = type,
            sel     = [convert_selection(sel1), convert_selection(sel2)],
            options = {'type' : type},
            result  = []
        )

    - You can use the convert_selection() function from the selection submodule to accept both strings \
        (that correspond to keys in the EMDA.selections dictionary) or AtomGroups.
    - Describe the adder by adding a docstring including at least the DESCRIPTION and USAGE sections, please.
"""


def add_distance(self, name, sel1, sel2, type="min"):
    """
    DESCRIPTION:
        This function outputs the minimum measured distance between the two input selections or coordinates or their combination.

    USAGE:
        EMDA.add_distance(name, sel1, sel2, type=['min' | 'cog' | 'com'])

    INPUT:
        - Name of the measurement
        - Two selections, which can contain more than one atom
        - type: min [default], com (center of mass) or cog (center of geometry)

    OUTPUT:
        - Shorter distance between sel1 and sel2 (in ang) or distances between COMs or COGs.


    """

    if type.lower() not in ("min", "max", "com", "cog"):
        raise NotAvailableOptionError

    # add the Measure dataclass to the measures list for the EMDA class
    self.measures[name] = self.Measure(
        name=name,
        type="distance",
        sel=[convert_selection(self, sel1), convert_selection(self, sel2)],
        options={"type": type},
        result=[],
    )

    return "Distance added!"


def add_angle(self, name, sel1, sel2, sel3, units="deg", domain=360):
    """
    DESCRIPTION:
        This functions measures the angle between 3 specified atoms and returns the value between 0 and 360 degrees.
        The input selections have to be single atoms.

    OPTIONS:
        - Name of the measurement
        - units: option for selecting the output units of the dihedral
            - degree
            - rad

        - domain: option for specifying the domain of the output measures
            - 180, pi: option for -180,180 domain
            - 360, 2pi: option for 0,360 domain. Default option

    INPUT:
        - Selection of four atoms in three different AtomGroups. They have to be input with the correct order

    OUTPUT:
        - Angle between the input atoms
    """

    for sel in (sel1, sel2, sel3):
        if isinstance(sel, AtomGroup):
            if len(sel) != 1:
                raise NotSingleAtomSelectionError
        elif isinstance(sel, str):
            if sel not in self.selections.keys():
                raise NotExistingSelectionError

    units = units.lower()
    domain = str(domain).lower()

    if units not in ("deg", "degree", "degrees", "rad", "radian", "radians"):
        units = "degree"

    if domain not in (180, "180", 360, "360", "pi", "2pi"):
        domain = "360"

    self.measures[name] = self.Measure(
        name=name,
        type="angle",
        sel=[
            convert_selection(self, sel1),
            convert_selection(self, sel2),
            convert_selection(self, sel3),
        ],
        options={"units": units, "domain": domain},
        result=[],
    )


def add_dihedral(self, name, sel1, sel2, sel3, sel4, units="degree", domain=360):
    """
    DESCRIPTION:
        This functions measures the dihedral angle between 4 specified atoms and returns the dihedral value between 0 and 360 degrees.
        The input selections have to be single atoms.

    OPTIONS:
        - units: option for selecting the output units of the dihedral
            - degree
            - rad

        - domain: option for specifying the domain of the output measures
            - 180, pi: option for -180,180 domain
            - 360, 2pi: option for 0,360 domain. Default option

    INPUT:
        - Name of the measurement
        - Selection of four atoms in four different AtomGroups. They have to be input with the correct order

    OUTPUT:
        - Dihedral angle between the input atoms
    """

    for sel in (sel1, sel2, sel3, sel4):
        if isinstance(sel, AtomGroup):
            if len(sel) != 1:
                raise NotSingleAtomSelectionError
        elif isinstance(sel, str):
            if sel not in self.selections.keys():
                raise NotExistingSelectionError

    units = units.lower()
    domain = str(domain).lower()

    if units not in ("deg", "degree", "degrees", "rad", "radian", "radians"):
        units = "degree"

    if domain not in (180, "180", 360, "360", "pi", "2pi"):
        domain = "360"

    self.measures[name] = self.Measure(
        name=name,
        type="dihedral",
        sel=[
            convert_selection(self, sel1),
            convert_selection(self, sel2),
            convert_selection(self, sel3),
            convert_selection(self, sel4),
        ],
        options={"units": units, "domain": domain},
        result=[],
    )


def add_planar_angle(self, name, sel1, sel2, units="deg", domain=360):
    """
    DESCRIPTION:
        This function measures the angle between two planes specified by three atoms each one and returns the angle.
        The input selections have to contain three atoms.

    OPTIONS:
        - Name of the measurement
        - units: option for selecting the output units of the dihedral
            - degree
            - rad

        - domain: option for specifying the domain of the output measures
            - 180, pi: option for -180,180 domain
            - 360, 2pi: option for 0,360 domain. Default option

    INPUT:
        - Selection of two sets of three atoms in two different AtomGroups.

    OUTPUT:
        - Angle between the input atoms
    """

    for sel in (sel1, sel2):
        if isinstance(sel, AtomGroup):
            if len(sel) != 3:
                raise NotSingleAtomSelectionError
        elif isinstance(sel, str):
            if sel not in self.selections.keys():
                raise NotExistingSelectionError
            else:
                if len(convert_selection(self, sel)) != 3:
                    raise NotSingleAtomSelectionError

    units = units.lower()
    domain = str(domain).lower()

    if units not in ("deg", "degree", "degrees", "rad", "radian", "radians"):
        units = "degree"

    if domain not in (180, "180", 360, "360", "pi", "2pi"):
        domain = "360"

    self.measures[name] = self.Measure(
        name=name,
        type="planar_angle",
        sel=[convert_selection(self, sel1), convert_selection(self, sel2)],
        options={"units": units, "domain": domain},
        result=[],
    )


def add_contacts(
    self,
    name,
    sel,
    sel_env=3,
    interactions="all",
    include_WAT=False,
    out_format="new",
    measure_distances=True,
):
    """
    DESCRIPTION:
        This function takes a Universe, a selection and a radius and returns the list of residues nearer than the specified radius.

    INPUT:
        - Name of the measurement
        - sel          -> selection of central atoms. It has to be an AtomGroup (only option for old engine) or the 'protein' string for
                            measuring distances between all contacting residues
        - sel_env      -> radius (in ang)
        - interactions -> type of interactions to be considered (all, polar, nonpolar, donorHbond, none). Custom
                            interactions can be also analysed by passing a list of residues names
        - out_format   -> [ '0.3'/'new'/'n' | '0.2'/'old'/'o' ] Format of the output. 'old' corresponds to the old logics (versions 0.0 to 0.2)
                            and is kept for compatibility reasons.

    OUTPUT:
        - List of dictionaries containing the name and number of all interacting residues
    """

    if isinstance(interactions, str):
        if interactions not in ("all", "polar", "nonpolar", "donorHbond", "none"):
            raise NotExistingInteractionError

        else:
            if interactions == "all":
                interactions = [
                    "ARG",
                    "HIS",
                    "HID",
                    "HIE",
                    "HIP",
                    "LYS",
                    "ASP",
                    "ASH",
                    "GLU",
                    "GLH",
                    "SER",
                    "THR",
                    "ASN",
                    "GLN",
                    "CYS",
                    "SEC",
                    "GLY",
                    "PRO",
                    "ALA",
                    "ILE",
                    "LEU",
                    "MET",
                    "PHE",
                    "TRP",
                    "TYR",
                    "VAL",
                ]

            elif interactions == "polar":
                interactions = [
                    "ARG",
                    "HIS",
                    "HID",
                    "HIE",
                    "HIP",
                    "LYS",
                    "ASP",
                    "ASH",
                    "GLU",
                    "GLH",
                    "SER",
                    "THR",
                    "ASN",
                    "GLN",
                    "CYS",
                    "SEC",
                    "TYR",
                ]

            elif interactions == "nonpolar":
                interactions = [
                    "CYS",
                    "SEC",
                    "GLY",
                    "PRO",
                    "ALA",
                    "ILE",
                    "LEU",
                    "MET",
                    "PHE",
                    "TRP",
                    "TYR",
                    "VAL",
                ]

            elif interactions == "donorHbond":
                interactions = [
                    "ARG",
                    "HID",
                    "HIE",
                    "HIP",
                    "LYS",
                    "ASH",
                    "GLH",
                    "SER",
                    "THR",
                    "ASN",
                    "GLN",
                    "CYS",
                    "SEC",
                    "GLY",
                    "PRO",
                    "TYR",
                ]

            elif interactions == "none":
                interactions = []

    elif isinstance(interactions, list):

        pass

    else:
        raise NotExistingInteractionError

    if include_WAT == True:
        interactions += ["WAT", "HOH"]

    if sel == "protein":
        mode = "protein"
        self.select("protein", sel, sel_type=None)
        sel = "protein"
        out_format = "new"

    elif isinstance(sel, AtomGroup) or sel in self.selections.keys():
        mode = "selection"

        if isinstance(sel, str):
            sel = convert_selection(self, sel)

        sel_env = self.universe.select_atoms(
            f"around {sel_env} group select", select=sel, updating=True
        )

        if str(out_format).lower() in ["0.2", "old", "o"]:
            out_format = "old"
            if measure_distances:
                print(
                    "Distances can not been calculated using the old output format. \
                        'measure_distances' has been set to False."
                )
                out_format = "new"

        elif str(out_format).lower() in ["0.3", "new", "n"]:
            out_format = "new"

        else:
            print(
                "The selected out format does not exist. The new format has been selected instead."
            )
            out_format = "new"

    self.measures[name] = self.Measure(
        name=name,
        type="contacts",
        sel=[convert_selection(self, sel), sel_env],
        options={
            "mode": mode,
            "interactions": interactions,
            "measure_dists": measure_distances,
            "out_format": out_format,
        },
        result=[],
    )


def add_RMSD(self, name, sel, ref=None, superposition=True):
    """
    DESCRIPTION:
        This function outputs the RMSD of a selection

    INPUT:
        - Name of the measurement
        - Selection
        - ref: selection of the reference universe. If not provided, the first frame will be used as the reference.
        - superposition [bool]: compute the RMSD of aligned

    OUTPUT:
        - Array of RMSDs of each frame against a reference
    """

    sel = convert_selection(self, sel)

    if isinstance(ref, type(None)):
        self.universe.trajectory[0]
        ref = sel.positions - sel.center_of_mass()

    elif isinstance(ref, AtomGroup):
        ref = ref.positions - ref.center_of_mass()

    self.measures[name] = self.Measure(
        name=name,
        type="RMSD",
        sel=[sel],
        options={"superposition": superposition, "ref": ref},
        result=[],
    )


def add_distWATbridge(self, name, sel1, sel2, sel1_rad=3, sel2_rad=3):
    """
    DESCRIPTION
        This function takes a Universe, two selections and the size of their environments and returns the nearest bridging water between the two selections and the distance to both of them.

    INPUT:
        - Name of the measurement
        - u            -> MDAnalysis Universe
        - sel1         -> selection of first set of central atoms. It has to be an AtomGroup
        - sel2         -> selection of second set of central atoms. It has to be an AtomGroup
        - sel1_rad      -> radius around the first set of central atoms (in ang)
        - sel2_rad      -> radius around the first set of central atoms (in ang)

    OUTPUT:
        - List of dictionaries containing the number of the bridging water and the smallest distance to each of the selection sets.
    """

    # sel1_rad, sel2_rad = sel1_env, sel2_env

    sel1_env = self.universe.select_atoms(
        "resname WAT and around %s group select" % sel1_rad,
        select=sel1,
        updating=True,
    )

    sel2_env = self.universe.select_atoms(
        "resname WAT and around %s group select" % sel2_rad,
        select=sel2,
        updating=True,
    )

    self.measures[name] = self.Measure(
        name=name,
        type="distWATbridge",
        sel=[
            convert_selection(self, sel1),
            convert_selection(self, sel2),
            convert_selection(self, sel1_env),
            convert_selection(self, sel2_env),
            convert_selection(self, sel1_rad),
            convert_selection(self, sel2_rad),
        ],
        options={},
        result=[],
    )


def add_pKa(
    self,
    name,
    excluded_ions=["Na+", "Cl-"],
    pka_ref="neutral",
    pdb_folder=".pka",
    keep_pdb=False,
    keep_pka=False,
):
    """
    DESCRIPTION:
        This function allows the prediction of the pKa using PROpKa3 of the protein for each frame.

    INPUT:
        - name:             name of the measurement
        - excluded_ions:    list of solvent ions names that belong to solvent. Default are Na+ and Cl-.
        - pka_ref:          reference to calculate pKa. Default is neutral. [ neutral | low-pH]
        - keep_pdb:         trigger for keeping generated pdbs. Default is False. [ True | False ]
        - keep_pka:         trigger for keeping generated .pka file. Default is False. [ True | False ]

    OUPUT:
        - Per-frame array of dicts with shape { residue : pKa }
    """

    excluded_ions = " or resname ".join(excluded_ions)

    sel = self.universe.select_atoms(
        "not (resname WAT or resname HOH or resname " + excluded_ions + ")"
    )

    self.measures[name] = self.Measure(
        name=name,
        type="pka",
        sel=[convert_selection(self, sel)],
        options={
            "pka_ref": pka_ref,
            "pdb_folder": pdb_folder,
            "keep_pdb": keep_pdb,
            "keep_pka": keep_pka,
        },
        result=[],
    )
