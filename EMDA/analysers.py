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
        min_val = val1-val2
        max_val = val1+val2
    
    elif mode.lower() in ('thres', 'threshold'):
        if val1 > val2:
            min_val, max_val =  val2, val1
        elif val1 < val2:
            min_val, max_val =  val1, val2

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


