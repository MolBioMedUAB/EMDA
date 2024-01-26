from .exceptions import NotCompatibleMeasureForAnalysisError,NotAvailableOptionError
from tqdm.autonotebook import tqdm

#from dataclasses import dataclass


def analyse_value(self, name, measure, val1, val2=0, mode='thres'):
    """
    DESCRIPTION:
        Analyser for checking if a value in the frame is between to given values. Threshold (a upper and lower (default is 0) limits) \
            and tolerance (val1 is the central value and val2 the +- range). If the value for each frame is within the set limits, \
            True will be returned in the corresponding step.

    OPTIONS:
        - val1, val2:   upper and lower limits (if threshold mode) or reference and tolerance values (if tolerance mode)
        - mode:         [ 'thres' | 'tol' ] Mode to create the satisfying range. Thres(hold) mode will return True if the value is \
                            between the given values, while tol(erance) mode will return True if the value is between val1+val2 and val1-val2.
    """

    if self.measures[measure].type not in ('distance', 'angle', 'dihedral', 'planar_angle'):
        raise NotCompatibleMeasureForAnalysisError
    

    if mode.lower() in ('tol', 'tolerance'):
        min_val, max_val = val1-val2, val1+val2
    
    elif mode.lower() in ('thres', 'threshold'):
        min_val, max_val = min(val1, val2), max(val1, val2)

    else :
        raise NotAvailableOptionError
    

    self.analyses[name] = self.Analysis(
        name = name,
        type = 'value',
        measure_name = measure,
        result = []
    )

    for result in self.measures[measure].result:
        self.analyses[name].result.append(
            min_val < result < max_val
        )



def contacts_frequency(self, name, percentage=False):
    """
    DESCRIPTION:
        
    """


    # identify if selection or whole protein contacts
    if isinstance(list(self.measurements.results[name][0].values())[0], float):
        contacts_type = 'sel'

    elif isinstance(list(self.measurements.results[name][0].values())[0], dict):
        contacts_type = 'prot'


    if contacts_type == 'prot':
        #results = []

        # create dict containing the residue name as key and a list as value. In this list, each contact in each frame will be stored
        total_contacts = {resid : [] for resid in list(self.measurements.results[name][0].keys())}

        for frame in self.measurements.results[name]:
            for resid in list(total_contacts.keys()):
                total_contacts[resid] += list(frame[resid].keys())

        contacts_freq = {}
        for residue in list(total_contacts.keys()):

            if percentage:
                contacts_freq[residue] = { residue_from_tot : total_contacts[residue].count(residue_from_tot)*100/len(self.measurements.results[name]) for residue_from_tot in list(set(list(total_contacts[residue])))}

            elif not percentage:
                contacts_freq[residue] = { residue_from_tot : total_contacts[residue].count(residue_from_tot) for residue_from_tot in list(set(list(total_contacts[residue])))}


    return contacts_freq
