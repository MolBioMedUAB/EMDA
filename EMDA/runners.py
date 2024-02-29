from .calculators import *


def run_distance(self, Measure, variant, replica):
    """
    DESCRIPTION:

    """

    Measure.result[variant][replica].append(
        calc_distance(
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[0]]),
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[1]]),
            Measure.options["type"])
    )


def run_angle(self, Measure, variant, replica):
    """
    DESCRIPTION:

    """

    Measure.result[variant][replica].append(
        calc_angle(
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[0]]),
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[1]]),
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[2]]),
            Measure.options["units"],
            Measure.options["domain"],
        )
    )


def run_dihedral(self, Measure, variant, replica):
    """
    DESCRIPTION:

    """

    Measure.result[variant][replica].append(
        calc_dihedral(
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[0]]),
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[1]]),
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[2]]),
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[3]]),
            Measure.options["units"],
            Measure.options["domain"],
        )
    )


def run_planar_angle(self, Measure, variant, replica):
    """
    DESCRIPTION:

    """

    Measure.result[variant][replica].append(
        calc_planar_angle(
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[0]]),
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[1]]),
            Measure.options["units"],
            Measure.options["domain"],
        )
    )


def run_pka(self, Measure, variant, replica):
    """
    DESCRIPTION:

    """

    Measure.result[variant][replica].append(
        calc_pka(
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[0]]),
            Measure.options["pka_ref"],
            Measure.options["pdb_folder"],
            Measure.options["keep_pdb"],
            Measure.options["keep_pka"],
        )
    )


def run_contacts(self, Measure, variant, replica):
    """
    DESCRIPTION:

    """

    Measure.result[variant][replica].append(
        calc_contacts_selection(
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[0]]),
            self.universe[variant][replica].select_atoms(
                f"around {Measure.sel[1]} group select",
                select=self.universe[variant][replica].select_atoms(self.selections[Measure.sel[0]])
            ),
            Measure.options["interactions"],
            Measure.options["measure_dists"],
            )
        )


def run_protein_contacts(self, Measure, variant, replica):
    """
    DESCRIPTION:

    """

    Measure.result[variant][replica].append(
        calc_contacts_protein(
            self.universe[variant][replica].select_atoms(Measure.sel[0]),
            Measure.sel[1],
            #self.universe[variant][replica].select_atoms(
            #    f"around {Measure.sel[1]} group select",
            #    select=self.universe[variant][replica].select_atoms(self.selections[Measure.sel[0]])
            #),
            measure_distances=Measure.options["measure_dists"],
            include_WAT=Measure.options["include_WAT"]
        )
    )


def run_RMSD(self, Measure, variant, replica):
    """
    DESCRIPTION:

    """

    Measure.result[variant][replica].append(
        calc_RMSD( 
            self.universe[variant][replica].select_atoms(self.selections[Measure.sel[0]]),
            Measure.options["ref"][variant][replica],
            Measure.options["superposition"]
        )
    )


def run_distWATbridge(self, Measure, variant, replica):
    """
    DESCRIPTION:

    """

    Measure.result[variant][replica].append(
        calc_distWATbridge(
            self.universe[variant][replica],
            self.universe[variant][replica].select_atoms(Measure.sel[0]),
            self.universe[variant][replica].select_atoms(Measure.sel[1]),
            self.universe[variant][replica].select_atoms(Measure.sel[2]),
            self.universe[variant][replica].select_atoms(Measure.sel[3]),
        )
    )
