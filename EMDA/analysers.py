from .exceptions import NotCompatibleMeasureForAnalysisError, NotAvailableOptionError, NotCompatibleContactsFormatError
from tqdm.autonotebook import tqdm

#from dataclasses import dataclass


def analyse_value(self, name, measure, val1, val2=0, mode='thres'):
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



def analyse_contacts_frequency(self, name, measure, percentage=False):
    """
    DESCRIPTION:
        Analyser for calculating the frequency (in absolute value or %) of the calculated contacts.

        This is a type of analysis that returns a dictiona

    OUTPUT:
        If the contacts is of the whole protein (so mode is protein), a dictionary in the result attribute of an Analysis class that contains each residue as key and 
            a dictionary with the residue to which contacts (as key) and the number of times it does or the % of frames (as value) as value.
        If the contacts is of a selection (so mode is selection), a dictionary containing the residue to which it contacts as key and the number of contacts or the 
            percentage as value.

    """

    if self.measures[measure].type not in ('contacts'):
        raise NotCompatibleMeasureForAnalysisError
    
    if self.measures[measure].options['out_format'] not in ('new'):
        raise NotCompatibleContactsFormatError

    if self.measures[measure].options['mode'] == 'protein':
        # create dict containing the residue name as key and a list as value. In this list, each contact in each frame will be stored
        total_contacts = {resid : [] for resid in list(self.measures[measure].result[0].keys())}

        for frame in self.measures[measure].result:
            for resid in list(total_contacts.keys()):
                total_contacts[resid] += list(frame[resid].keys())

        contacts_freq = {}
        for residue in list(total_contacts.keys()):

            if percentage:
                contacts_freq[residue] = { residue_from_tot : total_contacts[residue].count(residue_from_tot)*100/len(self.measures[measure].result) for residue_from_tot in list(set(list(total_contacts[residue])))}

            elif not percentage:
                contacts_freq[residue] = { residue_from_tot : total_contacts[residue].count(residue_from_tot) for residue_from_tot in list(set(list(total_contacts[residue])))}

    elif self.measures[measure].options['mode'] == 'selection':
        
        contacts_freq = {}
        for frame in self.measures[measure].result:
            for resid in list(frame.keys()):
                if resid in contacts_freq.keys():
                    contacts_freq[resid] += 1
                elif resid not in contacts_freq.keys():
                    contacts_freq[resid] = 1


    self.analyses[name] = self.Analysis(
        name = name,
        type = 'contacts_frequency',
        mode = self.measures[measure].options['mode'],
        measure_name = measure,
        result = contacts_freq
    )


#def analyse_NACs(self, name, analyses : list, percentage):