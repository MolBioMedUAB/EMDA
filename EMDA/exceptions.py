class NotSingleAtomSelectionError(Exception):
    """
    Raised when the input selection is not a single atom.
    """

    def __init__(self):
        Exception.__init__(self, "Only one atom has to be selected.")

    pass


class NotThreeAtomsSelectionError(Exception):
    """
    Raised when the input selection does not contain three atoms.
    """

    def __init__(self):
        Exception.__init__(self, "Three atoms have to be selected.")

    pass


class NotEqualListsLenghtError(Exception):
    """
    Raised when the 'sel_input' and the 'sel_types' list are not equal.
    """

    def __init__(self):
        Exception.__init__(self, "'sel_input' and the 'sel_types' list are not equal")

    pass


class NotEnoughAtomsSetectedError(Exception):
    """
    Raised when more than n atoms are required but not input.
    """

    def __init__(self):
        Exception.__init__(self, "Number of selected atoms is not enough")

    pass


class NotExistingMetalError(Exception):
    """
    Raised when the selected metal has no basis set (any of all the described) described.
    """

    def __init__(self):
        Exception.__init__(self, "The metal has no basis set available")

    pass


class NotExistingInteraction(Exception):
    """
    Raised when the selected interactions do not exist or are not available.
    """

    def __init__(self):
        Exception.__init__(self, "This type of interaction is not described. Available interactions are: 'all', 'polar', 'nonpolar', 'donorHbond' and 'none'. You can also add a custom list by inputing a list of residue names (in the three-letters coding)")

    pass


class OutputFormatNotAvailable(Exception):
    """
    Raised when the format of the output file is not available.
    """

    def __init__(self):
        Exception.__init__(self, "The output format is not available. Use JSON (.json or .jsn) or YAML (.yaml or .yml) instead.")

    pass


class NotAvailableOptionError(Exception):
    """
    Raised when an input option is not available.
    """

    def __init__(self):
        Exception.__init__(self, "One of the input options is not available. Revise the documentation of the function")

    pass

class EmptyMeasurementsError(Exception):
    """
    Raised when no measurment has been added when executing m.run_measure()
    """

    def __init__(self):
        Exception.__init__(self, "Add at least one measurement and run again.")

    pass
