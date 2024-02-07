# from ..exceptions import NotEqualListsLenghtError
from MDAnalysis.core.groups import AtomGroup


def selection(
    u, sel_input, sel_type=None, no_backbone=False, return_atomic_sel_string=False
):
    """
    DESCRIPTION
        This function takes an input number or name and type of selection (atom or residue or none) and transforms it into an MDAnalysis selection. The input can also be a list of numbers or name and it can combine both types of input by using a list of selection types (of the same lenght).

    INPUT
        - u: MDAnalysis' universe
        - sel_input: int, str, list or tuple of int, or list or tuple of str
        - sel_type: type of selection input as sel_input
            at_num    -> atom number
            at_name   -> atom name
            res_num   -> residue number
            res_name  -> residue name
            none/pipe -> pipes the selection command directly
        - no_backbone: bool. False is the default. True removes backbone atoms from selection
        - return_atomic_sel_string: Returns the selection string defined by atom

    OUTPUT
        - MDAnalysis selection as AtomGroup
    """

    if sel_type == None or sel_type.lower() in ("none", "pipe"):
        sel_string = sel_input

    elif isinstance(sel_input, list):
        if isinstance(sel_type, str):
            if sel_type.lower() == "at_num":
                sel_string = "bynum " + " or bynum ".join(str(at) for at in sel_input)
            elif sel_type.lower() == "at_name":
                sel_string = "name " + " or name ".join(sel_input)
            elif sel_type.lower() == "res_num":
                sel_string = "resid " + " or resid ".join(str(at) for at in sel_input)
            elif sel_type.lower() == "res_name":
                sel_string = "resname " + " or resname ".join(sel_input)
            else:
                sel_string = " or ".join(sel_input)

        elif isinstance(sel_type, list):
            if len(sel_input) != len(sel_type):
                pass  # raise NotEqualListsLenghtError

            else:
                sel_string = ""
                for inp, typ in zip(sel_input, sel_type):
                    if sel_string == "":
                        pass
                    else:
                        sel_string = "".join([sel_string, " or"])

                    if typ.lower() == "at_num":
                        sel_string = " ".join([sel_string, "bynum", str(inp)])
                    elif typ.lower() == "at_name":
                        sel_string = " ".join([sel_string, "name", str(inp)])
                    elif typ.lower() == "res_num":
                        sel_string = " ".join([sel_string, "resid", str(inp)])
                    elif typ.lower() == "res_name":
                        sel_string = " ".join([sel_string, "resname", str(inp)])
                    else:
                        sel_string = " ".join([sel_string, str(inp)])

    elif isinstance(sel_input, (str, int)):
        if sel_type.lower() == "at_num":
            sel_string = " ".join(["bynum", str(sel_input)])
        elif sel_type.lower() == "at_name":
            sel_string = " ".join(["name", str(sel_input)])
        elif sel_type.lower() == "res_num":
            sel_string = " ".join(["resid", str(sel_input)])
        elif sel_type.lower() == "res_name":
            sel_string = " ".join(["resname", str(sel_input)])

    if no_backbone == True:
        sel_string = (
            sel_string
            + " and not (name H or name N or name CA or name HA or name C or name O or name OXT or name H1 or name H2 or name H3)"
        )

    if return_atomic_sel_string == False:
        return u.select_atoms(sel_string)

    elif return_atomic_sel_string == True:
        return "index " + " or index ".join(list(u.select_atoms(sel_string).indices))


def convert_selection(self, sel):
    """
    DESCRIPTION
        Function for checking if the given selection is a string, so the AtomGroup has to be extracted from EMDA.selections, \
        or if it is an AtomGroup, so nothing has to be done.
    """

    if isinstance(sel, AtomGroup):
        return sel

    elif isinstance(sel, str):
        return self.selections[sel]
