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


class NotExistingInteractionError(Exception):
    """
    Raised when the selected interactions do not exist or are not available.
    """

    def __init__(self):
        Exception.__init__(
            self,
            "This type of interaction is not described. Available interactions are: 'all', 'polar', 'nonpolar', 'donorHbond' and 'none'. You can also add a custom list by inputing a list of residue names (in the three-letters coding)",
        )

    pass


class NotExistingSelectionError(Exception):
    """
    Raised when the selected selection do not exist or is not available.
    """

    def __init__(self):
        Exception.__init__(
            self,
            "This selection is not described. Create it and add the measure again.",
        )

    pass


class OutputFormatNotAvailableError(Exception):
    """
    Raised when the format of the output file is not available.
    """

    def __init__(self):
        Exception.__init__(
            self,
            "The output format is not available. Use JSON (.json or .jsn) or YAML (.yaml or .yml) instead.",
        )

    pass


class NotAvailableOptionError(Exception):
    """
    Raised when an input option is not available.
    """

    def __init__(self):
        Exception.__init__(
            self,
            "One of the input options is not available. Revise the documentation of the function",
        )

    pass


class EmptyMeasuresError(Exception):
    """
    Raised when no measurment has been added when executing m.run_measure()
    """

    def __init__(self):
        Exception.__init__(self, "Add at least one measurement and run again.")

    pass


class NotCompatibleMeasureForAnalysisError(Exception):
    """
    Raised when the measure to analyse is not of an accepted type by the analyser.
    """

    def __init__(self):
        Exception.__init__(
            self, "The input measure is not compatible with the chosen analysis."
        )

    pass


class NotCompatibleAnalysisForAnalysisError(Exception):
    """
    Raised when the analyses to analyse with a metaanalyser is not of an accepted type by the analyser.
    """

    def __init__(self):
        Exception.__init__(
            self, "The input analyses is not compatible with the chosen analysis."
        )

    pass


class NotCompatibleContactsFormatError(Exception):
    """
    Raised when the contacts format is not the new one.
    """

    def __init__(self):
        Exception.__init__(
            self,
            "The contacts format is not the new. Rerun the measure with the new format.",
        )

    pass


class NotEqualLenghtsError(Exception):
    """
    Raised when lists to compare do not have the same number of frames (so lists' lenghts).
    """

    def __init__(self, list_names, lenght):
        # list_names_str = ', '.join(list_names)
        Exception.__init__(
            self,
            f"{', '.join(list_names)} has/have not the same number of frames than the most common number ({lenght}).",
        )

    pass


class NotEnoughDataError(Exception):
    """
    Raised when the given data is not enough to perform the measure or analysis.
    """

    def __init__(self, minimum_data):
        Exception.__init__(
            self,
            f"At least {minimum_data} objects have to be given as input for performing the task.",
        )

    pass
