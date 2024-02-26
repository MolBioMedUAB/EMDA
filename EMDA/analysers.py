# Load package's functions
from .exceptions import (
    NotCompatibleMeasureForAnalysisError,
    NotAvailableOptionError,
    NotCompatibleContactsFormatError,
)
from .exceptions import (
    NotCompatibleAnalysisForAnalysisError,
    NotEqualLenghtsError,
    NotEnoughDataError,
)
from .tools import get_most_frequent

# from numpy import maximum as max

"""
TO BUILD in 0.2.0:
    - [X] NACs
    - [X] Contacts counter (how much contacts are in each frame)
    - [] Contacts checker --> checks if a contact is present in each frame
"""

"""
AVAILABLE:
    - analyse_value:
    - analyse_contacts_frequency:
    - analyse_contacts_amount:
    - analyse_NACs:
"""


def analyse_value(self, name, measure, val1, val2=0, mode="thres"):
    """
    DESCRIPTION:
        Analyser for checking if a value in the frame is between to given values. Threshold (a upper and lower (default is 0) limits) \
            and tolerance (val1 is the central value and val2 the +- range). If the value for each frame is within the set limits, \
            True will be returned in the corresponding step.

    OUTPUT:
        A frame-wise list containing boolean values depending on if the criteria have been satisfied or not.

    OPTIONS:
        - val1, val2:   upper and lower limits (if threshold mode) or reference and tolerance values (if tolerance mode)
        - mode:         [ 'thres' | 'tol' ] Mode to create the satisfying range. Thres(hold) mode will return True if the value is \
                            between the given values, while tol(erance) mode will return True if the value is between val1+val2 and val1-val2.
    """

    if self.measures[measure].type not in (
        "distance",
        "angle",
        "dihedral",
        "planar_angle",
    ):
        raise NotCompatibleMeasureForAnalysisError

    if mode.lower() in ("tol", "tolerance"):
        min_val, max_val = val1 - val2, val1 + val2

    elif mode.lower() in ("thres", "threshold"):
        min_val, max_val = min(val1, val2), max(val1, val2)

    else:
        raise NotAvailableOptionError

    results = []

    for result in self.measures[measure].result:
        results.append(min_val < result < max_val)

    self.analyses[name] = self.Analysis(
        name=name, type="value", measure_name=measure, result=results, options={}
    )


def analyse_contacts_frequency(
    self, name, measure, percentage=False, normalise_to_most_frequent=False
):
    """
    DESCRIPTION:
        Analyser for calculating the frequency (in absolute value or %) of the calculated contacts.

        This is a type of analysis that returns a dictiona

    OUTPUT:
        If the contacts is of the whole protein (so mode is protein), a dictionary in the result attribute of an Analysis class that contains each residue as key and
            a dictionary with the residue to which contacts (as key) and the number of times it does or the % of frames (as value) as value.
        If the contacts is of a selection (so mode is selection), a dictionary containing the residue to which it contacts as key and the number of contacts or the
            percentage as value.

    OPTIONS:
        - percentage: Returns the values in percentage
    """

    if self.measures[measure].type not in ("contacts"):
        raise NotCompatibleMeasureForAnalysisError

    if self.measures[measure].options["out_format"] not in ("new"):
        raise NotCompatibleContactsFormatError

    if self.measures[measure].options["mode"] == "protein":
        # create dict containing the residue name as key and a list as value. In this list, each contact in each frame will be stored
        total_contacts = {
            resid: [] for resid in list(self.measures[measure].result[0].keys())
        }

        for frame in self.measures[measure].result:
            for resid in list(total_contacts.keys()):
                total_contacts[resid] += list(frame[resid].keys())

        contacts_freq = {}
        for residue in list(total_contacts.keys()):

            if percentage:
                contacts_freq[residue] = {
                    residue_from_tot: total_contacts[residue].count(residue_from_tot)
                    * 100
                    / len(self.measures[measure].result)
                    for residue_from_tot in list(set(list(total_contacts[residue])))
                }

            elif not percentage:
                contacts_freq[residue] = {
                    residue_from_tot: total_contacts[residue].count(residue_from_tot)
                    for residue_from_tot in list(set(list(total_contacts[residue])))
                }

    elif self.measures[measure].options["mode"] == "selection":

        contacts_freq = {}
        for frame in self.measures[measure].result:
            for resid in list(frame.keys()):
                if resid in contacts_freq.keys():
                    contacts_freq[resid] += 1
                elif resid not in contacts_freq.keys():
                    contacts_freq[resid] = 1

        if percentage and not normalise_to_most_frequent:
            for residue in list(contacts_freq.keys()):
                contacts_freq[residue] = (
                    contacts_freq[residue] * 100 / len(self.measures[measure].result)
                )

        elif not percentage and normalise_to_most_frequent:
            for residue in list(contacts_freq.keys()):
                contacts_freq[residue] = contacts_freq[residue] / max(
                    contacts_freq.values()
                )

        elif percentage and normalise_to_most_frequent:
            for residue in list(contacts_freq.keys()):
                contacts_freq[residue] = (
                    contacts_freq[residue] * 100 / max(contacts_freq.values())
                )

    self.analyses[name] = self.Analysis(
        name=name,
        type="contacts_frequency",
        measure_name=measure,
        result=contacts_freq,
        options={"mode": self.measures[measure].options["mode"], "percentage": True},
    )


def analyse_contacts_amount(self, name, measure):
    """
    DESCRIPTION:
        Analyser for calculating how much contacts take place in each frame

    OUTPUT:
        If the contacts mode is protein:
            A frame-wise list containing a dictionary per frame containing the residue as key and the number of contacts of the residue as value

        If the contacts mode is selection:
            A frame-wise list containing the number of contacts for each frame.
    """

    if self.measures[measure].type not in ("contacts"):
        raise NotCompatibleMeasureForAnalysisError

    if self.measures[measure].options["mode"] == "protein":
        # create dict containing the residue name as key and a list as value. In this list, each contact in each frame will be stored
        contacts_amount = []

        for frame in self.measures[measure].result:
            contacts_amount.append(
                {resid: [] for resid in list(self.measures[measure].result[0].keys())}
            )
            for resid in list(contacts_amount[-1].keys()):
                contacts_amount[-1][resid] = len(frame[resid])

    elif self.measures[measure].options["mode"] == "selection":

        contacts_amount = []
        for frame in self.measures[measure].result:
            contacts_amount.append(len(frame))

    self.analyses[name] = self.Analysis(
        name=name,
        type="contacts_amount",
        options={
            "mode": self.measures[measure].options["mode"],
        },
        measure_name=measure,
        result=contacts_amount,
    )


def analyse_NACs(self, name, analyses: list, inverse: list = False):
    """
    DESCRIPTION:
        Metaanalyser (analyses two or more analyses) for combining boolean-output Analysis. It reads the boolean value corresponding to each analysis and returns True if all are True.


    OPTIONS:
        - name:         Name of the analysis
        - analyses:     List of analyses' names to analyse
        - inverse:      List of analyses' names which will be treated in the opposite way, so True will be False and viceversa.
    """

    # Check if input analyses are of the proper type
    if len(analyses) < 2:
        raise NotEnoughDataError(2)

    for analysis in analyses:

        if not isinstance(self.analyses[analysis].result, list):
            raise NotCompatibleAnalysisForAnalysisError

        else:
            if not isinstance(self.analyses[analysis].result[0], bool):
                raise NotCompatibleAnalysisForAnalysisError

    # Check if all analyses have the same number of frames
    length = get_most_frequent(
        [len(self.analyses[analysis].result) for analysis in analyses]
    )

    not_equal = [
        analysis
        for analysis in analyses
        if len(self.analyses[analysis].result) != length
    ]

    if len(not_equal) != 0:
        raise NotEqualLenghtsError(list_names=not_equal, lenght=length)

    if inverse != False:
        for inverse_ in inverse:
            if inverse_ not in analyses:
                print(
                    f"{inverse_} is not in analyses, so it's value will not be inverted."
                )

    self.analyses[name] = self.Analysis(
        name=name, type="NACs", measure_name=analyses, result=[], options={}
    )

    for frame in range(length):
        result_ = True
        for analysis in analyses:
            if inverse != False:
                if analysis in inverse:
                    result_ = result_ and not self.analyses[analysis].result[frame]
                else:
                    result_ = result_ and self.analyses[analysis].result[frame]
            else:
                result_ = result_ and self.analyses[analysis].result[frame]

        self.analyses[name].result.append(result_)
