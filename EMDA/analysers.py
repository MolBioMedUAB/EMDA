# Load package's functions
from .exceptions import (
    NotCompatibleMeasureForAnalysisError,
    NotAvailableOptionError,
    NotAvailableAnalysisError,
    NotAvailableMeasureError
)
from .exceptions import (
    NotCompatibleAnalysisForAnalysisError,
    NotEqualLenghtsError,
    NotEnoughDataError,
)
from .tools import get_most_frequent, get_dictionary_structure

from typing import Literal

import math
from scipy.stats import gaussian_kde
import numpy as np
import random

"""
AVAILABLE ANALYSERS:
    - analyse_value:
    - analyse_contacts_frequency:
    - analyse_contacts_amount:
    - analyse_contacts_presence:
    - analyse_NACs:
    - analyse_probability_density: 
"""

analyse_value_types = Literal['thres', 'threshold', 'tol', 'tolerance']
def analyse_value(self, name, measure, val1, val2=0, mode : analyse_value_types = "thres"):
    """
    DESCRIPTION:
        Analyser for checking if a value in the frame is between to given values. Threshold (a upper and lower (default is 0) limits) \
            and tolerance (val1 is the central value and val2 the +- range). If the value for each frame is within the set limits, \
            True will be returned in the corresponding step.

    COMPATIBLE DATA:
        - measures: distance, angle, dihedral, planar_angle
        - analyses: contacts_amount

    OUTPUT:
        A frame-wise list containing boolean values depending on if the criteria have been satisfied or not.

    ARGUMENTS:
        - val1, val2:   upper and lower limits (if threshold mode) or reference and tolerance values (if tolerance mode)
        - mode:         [ 'thres' | 'tol' ] Mode to create the satisfying range. Thres(hold) mode will return True if the value is \
                            between the given values, while tol(erance) mode will return True if the value is between val1+val2 and val1-val2.
    """

    # If measure, check type and return the Measure as obj
    if measure in list(self.measures.keys()):
        if self.measures[measure].type not in ("distance", "angle", "dihedral", "planar_angle"):
            raise NotCompatibleMeasureForAnalysisError
        
        obj = self.measures[measure]
    
    # If analysis, check type and return the Analysis as obj
    elif measure in list(self.analyses.keys()):
        if self.analyses[measure].type in ("contacts_amounts") and self.analyses[measure].options["mode"] in ("contacts"):
            pass
        else :
            raise NotCompatibleAnalysisForAnalysisError
        
        obj = self.analyses[measure]
        
    # define max and min values depending on mode,
    if mode.lower() in ("tol", "tolerance"):
        min_val, max_val = val1 - val2, val1 + val2

    elif mode.lower() in ("thres", "threshold"):
        min_val, max_val = min(val1, val2), max(val1, val2)

    else:
        raise NotAvailableOptionError
    
    # Analyse results
    results = get_dictionary_structure(obj.result, [])
    for variant in list(obj.result.keys()):
        for replica in list(obj.result[variant].keys()):
            for result in obj.result[variant][replica]:
                results[variant][replica].append(min_val < result < max_val)


    # save result
    self.analyses[name] = self.Analysis(
        name=name, type="value", measure_name=measure, result=results, options={}
    )


def analyse_contacts_frequency(self, name, measure, percentage : bool = False, normalise_to_most_frequent : bool = False):
    """
    DESCRIPTION:
        Analyser for calculating the frequency (in absolute value or %) of the calculated contacts.

        This is a type of analysis that returns a dictiona

    OUTPUT:
        If the contacts is of the whole protein (so mode is protein), a dictionary in the result attribute of an Analysis class that contains each residue as key and
            a dictionary with the residue to which contacts (as key) and the number of times it does or the % of frames (as value) as value.
        If the contacts is of a selection (so mode is selection), a dictionary containing the residue to which it contacts as key and the number of contacts or the
            percentage as value.

    ARGUMENTS:
        - percentage: Returns the values in percentage
    """

    if self.measures[measure].type == "contacts":
        contacts_freqs = get_dictionary_structure(self.measures[measure].result, {})
        for variant in list(self.measures[measure].result.keys()):
            for replica in list(self.measures[measure].result[variant].keys()):

                for frame in self.measures[measure].result[variant][replica]:
                    for resid in list(frame.keys()):
                        if resid in contacts_freqs[variant][replica].keys():
                            contacts_freqs[variant][replica][resid] += 1
                        elif resid not in contacts_freqs[variant][replica].keys():
                            contacts_freqs[variant][replica][resid] = 1

                if percentage and not normalise_to_most_frequent:
                    for residue in list(contacts_freqs[variant][replica].keys()):
                        contacts_freqs[variant][replica][residue] = (
                            contacts_freqs[variant][replica][residue] * 100 / len(self.measures[measure].result[variant][replica])
                        )

                elif not percentage and normalise_to_most_frequent:
                    for residue in list(contacts_freqs[variant][replica].keys()):
                        contacts_freqs[variant][replica][residue] = (
                            contacts_freqs[variant][replica][residue] / max(contacts_freqs[variant][replica].values())
                        )

                elif percentage and normalise_to_most_frequent:
                    for residue in list(contacts_freqs[variant][replica].keys()):
                        contacts_freqs[variant][replica][residue] = (
                            contacts_freqs[variant][replica][residue] * 100/ max(contacts_freqs[variant][replica].values())
                        )
        
        
    elif self.measures[measure].type == "per_residue_contacts":
        contacts_freqs = get_dictionary_structure(self.measures[measure].result, {})
        for variant in list(self.measures[measure].result.keys()):
            for replica in list(self.measures[measure].result[variant].keys()):

                # create dict containing the residue name as key and a list as value. In this list, each contact in each frame will be stored
                total_contacts = {
                    resid: [] for resid in list(self.measures[measure].result[variant][replica][0].keys())
                }

                for frame in self.measures[measure].result[variant][replica]:
                    for resid in list(total_contacts.keys()):
                        total_contacts[resid] += list(frame[resid].keys())

                contacts_freqs[variant][replica] = {}
                for residue in list(total_contacts.keys()):
                
                    if percentage and not normalise_to_most_frequent:
                        contacts_freqs[variant][replica][residue] = {
                            residue_from_tot : total_contacts[residue].count(residue_from_tot) * 100 / len(self.measures[measure].result[variant][replica]) for residue_from_tot in list(set(list(total_contacts[residue])))
                        }

                    elif not percentage or normalise_to_most_frequent:
                        contacts_freqs[variant][replica][residue] = {
                            residue_from_tot: total_contacts[residue].count(residue_from_tot) for residue_from_tot in list(set(list(total_contacts[residue])))
                        }

                if normalise_to_most_frequent:
                    if percentage:
                        max_value = 0
                        for residue in list(contacts_freqs[variant][replica].keys()):
                            if max_value < max(list(contacts_freqs[variant][replica][residue].values())): 
                                max_value = max(list(contacts_freqs[variant][replica][residue].values()))

                        for residue in list(contacts_freqs[variant][replica].keys()):
                            for residue_ in list(contacts_freqs[variant][replica][residue].keys()):
                                contacts_freqs[variant][replica][residue][residue_] = contacts_freqs[variant][replica][residue][residue_] * 100 / max_value

                    elif not percentage:
                        max_value = 0
                        for residue in list(contacts_freqs[variant][replica].keys()):
                            if max_value < max(list(contacts_freqs[variant][replica][residue].values())): 
                                max_value = max(list(contacts_freqs[variant][replica][residue].values()))


                        #max(list(contacts_freqs[variant][replica].values()))
                        for residue in list(contacts_freqs[variant][replica].keys()):
                            for residue_ in list(contacts_freqs[variant][replica][residue].keys()):
                                contacts_freqs[variant][replica][residue][residue_] = contacts_freqs[variant][replica][residue][residue_] / max_value

                    #if percentage and not normalise_to_most_frequent:
                    #    for residue in list(contacts_freqs[variant][replica].keys()):
                    #        contacts_freqs[variant][replica][residue] = (
                    #            contacts_freqs[variant][replica][residue] * 100 / len(self.measures[measure].result[variant][replica])
                    #        )
#
                    #elif not percentage and normalise_to_most_frequent:
                    #    for residue in list(contacts_freqs[variant][replica].keys()):
                    #        contacts_freqs[variant][replica][residue] = (
                    #            contacts_freqs[variant][replica][residue] / max(contacts_freqs[variant][replica].values())
                    #        )
#
                    #elif percentage and normalise_to_most_frequent:
                    #    for residue in list(contacts_freqs[variant][replica].keys()):
                    #        contacts_freqs[variant][replica][residue] = (
                    #            contacts_freqs[variant][replica][residue] * 100/ max(contacts_freqs[variant][replica].values())
                    #        )
#
    else :#
        raise NotCompatibleMeasureForAnalysisError

    # per_residue_contacts-related code
    #if self.measures[measure].options["mode"] == "protein":
    #    # create dict containing the residue name as key and a list as value. In this list, each contact in each frame will be stored
    #    total_contacts = {
    #        resid: [] for resid in list(self.measures[measure].result[0].keys())
    #    }
#
    #    for frame in self.measures[measure].result:
    #        for resid in list(total_contacts.keys()):
    #            total_contacts[resid] += list(frame[resid].keys())
#
    #    contacts_freq = {}
    #    for residue in list(total_contacts.keys()):
#
    #        if percentage:
    #            contacts_freq[residue] = {
    #                residue_from_tot: total_contacts[residue].count(residue_from_tot)
    #                * 100
    #                / len(self.measures[measure].result)
    #                for residue_from_tot in list(set(list(total_contacts[residue])))
    #            }
#
    #        elif not percentage:
    #            contacts_freq[residue] = {
    #                residue_from_tot: total_contacts[residue].count(residue_from_tot)
    #                for residue_from_tot in list(set(list(total_contacts[residue])))
    #            }


    
        

    self.analyses[name] = self.Analysis(
        name=name,
        type="contacts_frequency",
        measure_name=measure,
        result=contacts_freqs,
        options={"mode": self.measures[measure].type, "percentage": percentage},
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

    if self.measures[measure].type  == "contacts":
        
        contacts_amount = get_dictionary_structure(self.measures[measure].result, [])
        for variant in list(self.measures[measure].result.keys()):
            for replica in list(self.measures[measure].result[variant].keys()):

                for frame in self.measures[measure].result[variant][replica]:
                    contacts_amount[variant][replica].append(len(frame))

    elif self.measures[measure].type == "per_residue_contacts":

        # create dict containing the residue name as key and a list as value insite the variant/replica dict. In this list, each contact in each frame will be stored
        contacts_amount = get_dictionary_structure(self.measures[measure].result, [])
        for variant in list(self.measures[measure].result.keys()):
            for replica in list(self.measures[measure].result[variant].keys()):

                for frame in self.measures[measure].result[variant][replica]:
                    contacts_amount[variant][replica].append(
                        {resid: [] for resid in list(self.measures[measure].result[variant][replica][0].keys())}
                    )
                    for resid in list(contacts_amount[variant][replica][-1].keys()):
                        contacts_amount[variant][replica][-1][resid] = len(frame[resid])

    else :
        raise NotCompatibleMeasureForAnalysisError


    self.analyses[name] = self.Analysis(
        name=name,
        type="contacts_amount",
        options={
            "mode": self.measures[measure].type,
        },
        measure_name=measure,
        result=contacts_amount,
    )

analyse_contacts_presence_mode_types = Literal['all', 'any']
def analyse_contacts_presence(self, name, measure, contact, mode : analyse_contacts_presence_mode_types = 'all'):
    """
    DESCRIPTION:
        Analyser for checking if a contacts or contacts is/are present in a contacts analysis. If the contact is present, a True is returned-.
        Currently only available for contact and not for protein_contact

    ARGUMENTS:
        - contact:  str or link containing the residue's number of the residue to be contacted
        - mode:     'all' returns
        requires all specified contacts to be present. 'any' requires at least one contact to be present.
    """

    # check the type of the measure type
    if self.measures[measure].type not in ('contacts'):
        raise NotCompatibleMeasureForAnalysisError
    
    # fix contact list if str is input
    if isinstance(contact, (str, int)):
        contact = [str(contact)]
    
    result = get_dictionary_structure(self.measures[measure].result, [])
    for variant in list(self.measures[measure].result.keys()):
        for replica in list(self.measures[measure].result[variant].keys()):
            for frame in range(len(self.measures[measure].result[variant][replica])):
                contacting = [ contacting_resid[3:] for contacting_resid in list(self.measures[measure].result[variant][replica][frame].keys())]
                result_ = []
                for c in contact:
                    if c in contacting:
                        result_.append(True)
                    else :
                        result_.append(False)

                if mode == 'all':
                    result[variant][replica].append(all(result_))
                elif mode == 'any':
                    result[variant][replica].append(any(result_))

            
    self.analyses[name] = self.Analysis(
        name=name,
        type="contacts_presence",
        options={
            "mode": self.measures[measure].type,
        },
        measure_name=measure,
        result=result,
    )


def analyse_NACs(self, name, analyses : list, merge_replicas : bool = False, invert : list = False):
    """
    DESCRIPTION:
        Metaanalyser (analyses two or more analyses) for combining boolean-output Analysis. It reads the boolean value corresponding to each analysis and returns True if all are True.


    ARGUMENTS:
        - name:         Name of the analysis
        - analyses:     List of analyses' names to analyse
        - inverse:      List of analyses' names which will be treated in the opposite way, so True will be False and viceversa.
    """

    # Check if input analyses are of the proper type
    if len(analyses) < 2:
        raise NotEnoughDataError(2) 

    for analysis in analyses:
        if self.analyses[analysis].type not in ("value", "contacts_presence"):
            raise NotCompatibleAnalysisForAnalysisError

    if invert != False:
        for invert_ in invert:
            if invert_ not in analyses:
                print(
                    f"{invert_} is not in analyses, so it's value will not be inverted."
                )

    if not merge_replicas:

        lengths = get_dictionary_structure(self.analyses[analyses[0]].result, {})

        for variant in list(self.analyses[analyses[0]].result.keys()):
            for replica in list(self.analyses[analyses[0]].result[variant].keys()):
                # Check if all analyses in same replica have the same number of frames
                lengths[variant][replica] = get_most_frequent([ len(self.analyses[analysis].result[variant][replica]) for analysis in analyses ])
                
                # Check those replicas with different length
                not_equal = [ analysis for analysis in analyses if len(self.analyses[analysis].result[variant][replica]) != lengths[variant][replica] ]


        if len(not_equal) != 0:
            raise NotEqualLenghtsError(list_names=not_equal)



        results = get_dictionary_structure(self.analyses[analyses[0]].result, []) 
        for variant in list(self.analyses[analyses[0]].result.keys()):
            for replica in list(self.analyses[analyses[0]].result[variant].keys()):
                result_ = True
                for result_ix in range(len(self.analyses[analyses[0]].result[variant][replica])):
                    for analysis in analyses:
                        # check if analysis name is false or different
                        if invert == False:
                            result_ = result_ and self.analyses[analysis].result[variant][replica][result_ix]
                        
                        else :
                            # check if analysis name is in invert or not
                            if analysis not in invert or not invert:
                                result_ = result_ and self.analyses[analysis].result[variant][replica][result_ix]
                            elif analysis in invert:
                                result_ = result_ and not self.analyses[analysis].result[variant][replica][result_ix]
                        
                    results[variant][replica].append(result_)
                    

    self.analyses[name] = self.Analysis(
        name=name,
        type="NACs",
        measure_name=analyses,
        result=results,
        options = {"merge_replicas" : merge_replicas}
    )


def analyse_probability_density(self, name, measures, bw_method = 'scott', get_basins : bool = True, num_of_points = None, print_results : bool = False):
    """
    DESCRIPTION:
        Analyser for getting the probability map for a certain event (distance, for instance).

    ARGUMENTS:
        - measures:     list of two measures to be taken into consideration
        - bw_method:    [ scott | silverman | float ] method to calculate the estimator bandwidth. More info at https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html
        - get_basins:   gets the minima get_density_v4.py script.

    SOURCE:
        Code from Bruno Victor's
    """

    if self.measures[measures[0]].type not in ('distance', 'angle', 'dihedral', 'planar_angle', 'RMSD'):#, 'contacts_amount'):
        raise NotCompatibleMeasureForAnalysisError(measure=measures[0])
    
    if self.measures[measures[1]].type not in ('distance', 'angle', 'dihedral', 'planar_angle', 'RMSD'):#, 'contacts_amount'):
        raise NotCompatibleMeasureForAnalysisError(measure=measures[1])


    # Find maximum value in a 2D array of n x n elements
    def calcPmax(PDF, n):
        pmax = 0
        for i in range(n):
            for j in range(n):
                if PDF[i][j] > pmax:
                    pmax = PDF[i][j]
        return pmax

    # Write the landscape to SDTOUT
    def writeLscape(n, PDF, xbin_width, ybin_width, pmax, X, Y, pdf_hist, x):
        lscape = []

        for i in range(n):
            for j in range(n):
                if PDF[i][j] > 0:
                    e = -math.log(PDF[i][j]/pmax)
                else:
                    e = 10000
                
                if print_results:
                    print(X[i][j], Y[i][j], PDF[i][j],
                        pdf_hist[i][j] / len(x) / (xbin_width * ybin_width),
                        e
                    )
                
                lscape.append([
                    X[i][j], Y[i][j], PDF[i][j],
                    pdf_hist[i][j] / len(x) / (xbin_width * ybin_width),
                    e
                ])
            if print_results:
                print()

        return lscape

    # Return xxs and yys axis from the X and Y numpy array that is used to
    # obtain the discretized PDF function
    def genAxes(X, Y):
        xxs = []
        yys = []
        for i in X:
            xxs.append(i[0])
        for i in Y[0]:
            yys.append(i)
        return xxs, yys


    # Getbasins - an index value will be atributed to each node in a n x n
    # grid
    def getBasins(PDF, xxs, yys, x, y, nData):
        if nData == len(x):
            ndx_sample = [i for i in range(len(x))]
        else:
            random.seed(23)
            ndx_sample = random.sample([i for i in range(len(x))], nData)
        all_min = []
        for i in ndx_sample:
            all_min.append(getSingleBasin(i, PDF, xxs, yys, x, y))
        return all_min, ndx_sample


    # Several getBasins (used by getBasins main function)
    def getSingleBasin(ndx, PDF, xxs, yys, x, y):
        basins = []
        first_node = findClosestNode(xxs, yys, x[ndx], y[ndx])
        basin = steepD(PDF, first_node)
        return basin

    # Run steepest descent
    def steepD(PDF, node):
        cnt = False
        x = node[0]
        y = node[1]
        new_node = [x, y]
        ##
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                try:
                    if PDF[x+i][y+j] > PDF[new_node[0]][new_node[1]]:
                        cnt = True
                        new_node = [x+i, y+j]
                except IndexError:
                    None
        ##
        if cnt is not False:
            return steepD(PDF, new_node)
        else:
            return new_node


    # Find closest bin in xxs and yys
    def findClosestNode(xxs, yys, x, y):
        i = findClosest1D(xxs, x)
        j = findClosest1D(yys, y)
        return [i, j]


    # Find closest bin in 1D
    def findClosest1D(arr, val, min_dist = 1000000):
        closest = 0
        for i in range(len(arr)):
            d = abs(arr[i] - val)
            if d < min_dist:
                min_dist = d
                closest = i
        return closest


    # Gen a list of minumum values without repetitions
    def genMinList(all_min):
        min_list = []
        for mini in all_min:
            if mini not in min_list:
                min_list.append(mini)
        min_probs = calcMinProbs(all_min, min_list)
        new_min_list = [mini for _, mini in sorted(zip(min_probs, min_list),
                                                reverse=True)]
        new_min_probs = sorted(min_probs, reverse=True)
        return new_min_list, new_min_probs


    # Calculate probabilities of each minimum
    def calcMinProbs(all_min, min_list):
        probs = []
        n = len(all_min)
        for i in min_list:
            cnt = 0
            for j in all_min:
                if i == j:
                    cnt += 1
            probs.append(cnt/n)
        return probs


    def run(measure1, measure2, bw_method = bw_method, get_basins = get_basins, num_of_points = None):
        """
        DESCRIPTION:
            Function for running the analysis for one variant and replica
        """
        
        if num_of_points == None:
            num_of_points = len(measure1)

        npoints = complex(num_of_points)

        try:
            bwm = float(bw_method)
        except ValueError:
            bwm = bw_method

        xmin = min(measure1)
        xmax = max(measure1)
        ymin = min(measure2)
        ymax = max(measure2)
        values = np.vstack([measure1, measure2])
        xbin_width = (xmax - xmin) / int(num_of_points)
        ybin_width = (ymax - ymin) / int(num_of_points)
        
        kernel = gaussian_kde(values, bw_method=bwm)
        # Save kernel variable into a binary file
        #pickle.dump(kernel, open("Pickle_dump-kernel.p", "wb" ) )

        X, Y = np.mgrid[xmin:xmax:npoints, ymin:ymax:npoints]
        positions = np.vstack([X.ravel(), Y.ravel()])
        PDF = np.reshape(kernel(positions).T, X.shape)

        edgesx = np.arange(xmin, xmax + xbin_width, xbin_width)
        edgesy = np.arange(ymin, ymax + ybin_width, ybin_width)
        pdf_hist, edgesx, edgesy = np.histogram2d(measure1, measure2, bins=(edgesx, edgesy))

        pmax = calcPmax(PDF, int(num_of_points))

        # OUTPUT landscape
        lscape = writeLscape(int(num_of_points), PDF, xbin_width, ybin_width, pmax, X=X, Y=Y, pdf_hist=pdf_hist, x=measure1)

        data = []
        mins = []
        if get_basins:
            xxs, yys = genAxes(X, Y)
            
            all_min, ndx_sample = getBasins(PDF, xxs, yys, measure1, measure2, int(num_of_points))
            
            min_list, min_probs = genMinList(all_min)

            # OUTPUT information for all data points
            for i in range(len(all_min)):
                for j in range(len(min_list)):
                    if all_min[i] == min_list[j]:
                        node = findClosestNode(xxs, yys, measure1[ndx_sample[i]],
                                            measure2[ndx_sample[i]])
                        
                        if print_results:
                            print("DATA", j, measure1[ndx_sample[i]], measure2[ndx_sample[i]],
                                -math.log(PDF[node[0]][node[1]]/pmax))
                        
                        data.append([j, measure1[ndx_sample[i]], measure2[ndx_sample[i]],
                            -math.log(PDF[node[0]][node[1]]/pmax)])

            # OUTPUT information for all basins
            for i in range(len(min_list)):

                if print_results:
                    print("MIN", i,
                        xxs[min_list[i][0]],
                        yys[min_list[i][1]], min_probs[i],
                        PDF[min_list[i][0]][min_list[i][1]],
                        -math.log(PDF[min_list[i][0]][min_list[i][1]]/pmax))
                
                mins.append([
                    i,
                    xxs[min_list[i][0]],
                    yys[min_list[i][1]], min_probs[i],
                    PDF[min_list[i][0]][min_list[i][1]],
                    -math.log(PDF[min_list[i][0]][min_list[i][1]]/pmax)
                ])
                
        return lscape, data, mins
        
    result_ = get_dictionary_structure(self.measures[measures[0]].result, {'lscape' : [], 'data' : [], 'mins' : []})
    #datas   = get_dictionary_structure(self.measures[measures[0]].result, [])
    #mins    = get_dictionary_structure(self.measures[measures[0]].result, [])
    for variant in list(self.measures[measures[0]].result.keys()):
        for replica in list(self.measures[measures[0]].result[variant].keys()):
            result_[variant][replica]['lscape'], result_[variant][replica]['data'], result_[variant][replica]['mins'] = run(
                measure1=self.measures[measures[0]].result[variant][replica],
                measure2=self.measures[measures[1]].result[variant][replica],
                bw_method='scott',
                get_basins=get_basins,
                num_of_points=num_of_points,
            )

    self.analyses[name] = self.Analysis(
        name=name,
        type="pdf",
        measure_name=measures,
        result=result_,
        options = {
            "bw_method"         : bw_method,
            "get_basins"        : get_basins,
            "num_of_points"     : num_of_points,
            "measure_types"     : [self.measures[measures[0]].type, self.measures[measures[1]].type],
            "selection_names"   : [','.join(self.measures[measures[0]].sel), ','.join(self.measures[measures[1]].sel)]
        }
    )