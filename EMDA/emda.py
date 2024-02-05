from dataclasses import dataclass, field, is_dataclass
import pickle

#from typing import Literal
#from typing import Union
#from types import NoneType

# load MDAnalysis' universe class
from MDAnalysis import Universe
#from MDAnalysis.core.groups import AtomGroup

# load internal EMDA classes and functions
from .selection import selection
from .adders import *
from .runners import *
from .analysers import *
from .plotters import *
from .tools import in_notebook

# load custom exceptions
from .exceptions import EmptyMeasuresError

from tqdm.autonotebook import tqdm

#from metaclass import add_adders
#class EMDA(metaclass=add_adders):

"""
IDEAS:
    - [] Make parameters and trajectory not only to strings but to dict with mutants and replicas.
"""


class EMDA:

    def __init__(self, parameters, trajectory):

        self.parameters = parameters
        self.trajectory = trajectory
        self.universe  = Universe(parameters, trajectory)
        print('Trajectory has been loaded!')
        self.selections = {}
        self.measures   = {}
        self.analyses   = {}

        
        # Automatically add all imported functions from adders.py and from analysers.py as EMDA methods
        external_functions  = [func for func in globals() if callable(globals()[func]) and func.startswith("add_")]
        external_functions += [func for func in globals() if callable(globals()[func]) and func.startswith("analyse_")]
        for func_name in external_functions:
            setattr(EMDA, func_name, globals()[func_name])

    @dataclass
    class Measure:
        name    : str
        type    : str
        sel     : list
        options : dict
        result  : list

        def __str__(self) -> str:
            if len(self.result) == 0:
                status = 'Not calculated'
            else :
                status = 'Calculated'

            print_  = f"Measure dataclass with:\n"
            print_ += f"\tName:   {self.name}\n"
            print_ += f"\tType:   {self.type}\n"
            print_ += f"\tSel:    {self.sel}\n"
            print_ += f"\tStatus: {status}\n"
        
            return print_
        
        def __repr__(self):
            return self.__str__()
        
        def plot(self):
            if self.type in ('distance', 'angle', 'dihedral', 'RMSD', 'planar_angle'):
                units = {
                    'distance' : '(Å)',
                    'angle' : '(º)',
                    'dihedral' : '(º)',
                    'planar_angle' : '(º)',
                    'RMSD' : '(Å)',
                         }

                plt.plot(self.result)
                plt.ylabel(' '.join(self.type.split('_')).capitalize() + ' ' + units[self.type])
                plt.xlabel('Frame')

                plt.show()
                plt.close()
            else :
                print(f"This method is still not available for {type} measures.")

        
    @dataclass
    class Analysis:
        name : str
        type : str
        measure_name : str
        result : list
        options : dict = field(default_factory={})
        #options : dict
        #mode : Union[str, NoneType] = None


        def __str__(self) -> str:
            if len(self.result) == 0:
                status = 'Not calculated'
            else :
                status = 'Calculated'

            print_  = f"Analysis dataclass with:\n"
            print_ += f"\tName:   {self.name}\n"
            print_ += f"\tType:   {self.type}\n"
            print_ += f"\tRelated mesure:    {self.measure_name}\n"
            print_ += f"\tStatus: {status}\n"

            return print_
        
        def __repr__(self):
            return self.__str__()

    # function for creating selections (AtomGroups) as a dictionary inside EMDA class
    def select(self, name, sel_input, sel_type=None, no_backbone=False, return_atomic_sel_string=False):
        self.selections[name] = selection(self.universe, sel_input, sel_type=sel_type, no_backbone=no_backbone, return_atomic_sel_string=return_atomic_sel_string)


    def print_available_adders(self):

        import inspect
        
        adder_methods = [method for method in dir(self) if method.startswith('add_')]
        
        for method_name in adder_methods:
            method = getattr(self, method_name)

            docstring = inspect.getdoc(method).split('\n')

            # Extract DESCRIPTION and USAGE sections from docstring
            description_section = ''
            usage_section       = ''
            for l in range(len(docstring)):
                if 'DESCRIPTION' in docstring[l]:
                    description_section += docstring[l]  + '\n'
                    l_ = l
                    while True:
                        l_ += 1
                        if docstring[l_] == '':
                            break
                        else :
                            description_section += docstring[l_]  + '\n'
                
                if 'USAGE' in docstring[l]:
                    usage_section += docstring[l] + '\n'
                    l_ = l
                    while True:
                        l_ += 1
                        if docstring[l_] == '':
                            break
                        else :
                            usage_section += docstring[l_]  + '\n'

            if description_section == '':
                description_section = 'No description available.\n'
            if usage_section == '':
                usage_section = 'No usage available.\n'

            print(f"\033[1mMethod: {method_name}\033[0m\nHelp:\n{description_section}\n{usage_section}\n")

        print("Use '>>> help(EMDA.add_***)' to get the complete information of an adder.")


    def run(self, exclude=None, run_only=None, recalculate=False, step=1, start=1, end=-1,):
        """ 
        DESCRIPTION:
            Run all the measurements configured in self.measures

        OPTIONS:
            - exclude:      skip measures with the given name. Ignored if used with run_only.
            - run_only:     run only the measures passed by list of names. If used with exclude, exclude will be ignored.
            - recalculate:  [True | False | (List | str)] True for resetting all precalculated measures, list or str with
                            the specific measures to recalculate.
            - step:         Frames to jump during the analysis. Default is 1, so all the trajectory will be analysed.
            - start:        First frame to start the analysis. Default is 0.
            - end:          Last frame to analyse (included). Default is last frame of trajectory.
        """
        
        # Check that there is at least one measure set
        if len(self.measures) == 0:
            raise EmptyMeasuresError
        
        # If last frame is -1 (the default), last in trajectory is chosen.
        if end == -1: end = len(self.universe.trajectory)

        # Convert exclude to list to append precalculated measures if recalculate is False.
        ## If it is True, set measure's result as empty list, so it is overwritten.
        ## If a measure or list of measures is given, their result list will be reset as [].
        if exclude == None: exclude = []

        if isinstance(recalculate, str):
                recalculate = [recalculate]

        for measure in list(self.measures.keys()):
            if len(self.measures[measure].result) > 0:
                if isinstance(recalculate, bool):
                    if recalculate:
                        self.measures[measure].result = []
                    elif not recalculate:
                        exclude += [measure]

                elif isinstance(recalculate, list):
                    if measure in recalculate:
                        self.measures[measure].result = []


        # Check run_only. If run_only is used, the measure set will be the run_only list. Conversely, measures will be set as 
        ## all measures except exclude (if none, it is converted to empty list in previous codeblock)
        if run_only != None:
            if isinstance(run_only, str):
                measures = set([run_only])
            elif isinstance(run_only, list):
                measures = set(run_only)
        else :
            measures = (set(self.measures.keys()) - set(exclude))

        # trajectory cycle
        first_cycle = True
        for ts in tqdm(self.universe.trajectory[start-1:end:step], desc='Measuring', unit='Frame'):
            
            # measures cycle
            for measure in measures:
                """
                TO-DO:
                    Look for an automatic way to run the proper runner depending on the type set in the Measure. \
                    All runners should be named like run_${type_of_measure}, so maybe there is a chance to automate it.
                """

                if self.measures[measure].type == 'distance':
                    run_distance(self.measures[measure])

                elif self.measures[measure].type == 'angle':
                    run_angle(self.measures[measure])

                elif self.measures[measure].type == 'dihedral':
                    run_dihedral(self.measures[measure])

                elif self.measures[measure].type == 'planar_angle':
                    run_planar_angle(self.measures[measure])

                elif self.measures[measure].type == 'contacts':
                    run_contacts(self.measures[measure])

                elif self.measures[measure].type == 'RMSD':
                    run_RMSD(self.measures[measure])
                
                elif self.measures[measure].type == 'distWATbridge':
                    run_distWATbridge(self.measures[measure])

                elif self.measures[measure].type == 'pka':
                    if first_cycle:
                        check_folder(self.measures[measure].options['pdb_folder'])

                    run_pka(self.measures[measure])

                else :
                    if first_cycle:
                        print(f"The {measure} type is not available. If you need, you can create it.")

            if first_cycle:
                first_cycle = False


    def save_result(self, name, out_name=None, ):
        
        if out_name == None: out_name = name + '.pickle'
        
        if name in list(self.analyses.keys()):
            with open('.'.joint(out_name.split('.')[:-1]) + '_analysis' + out_name.split('.')[-1], 'wb') as handle:
                pickle.dump(self.analyses[name].result, handle, protocol=2)

            print(f"{name} analysis' result has been saved as {out_name}!")

        elif name in list(self.measures.keys()):
            with open('.'.joint(out_name.split('.')[:-1]) + '_measure' + out_name.split('.')[-1], 'wb') as handle:
                pickle.dump(self.measures[name].result, handle, protocol=2)

            print(f"{name} analysis' result has been saved as {out_name}!")

        else :
            raise KeyError(f"{name} is not an available measure nor analysis.")


    def read_result(self, filename, name,  type = "dataclass"):

        if is_dataclass(name) and type.lower() in ["d", "dataclass"]:
            with open(filename, 'rb') as handle:
                name.result = pickle.load(handle)

        else :
            try :
                if type.lower() in ["m", "measure"]:
                    with open(filename, 'rb') as handle:
                        self.measures[name].result = pickle.load(handle)

                if type.lower() in ["a", "analysis"]:
                    with open(filename, 'rb') as handle:
                        self.analyses[name].result = pickle.load(handle)

            except KeyError:
                raise KeyError(f"{name} is not an available measure nor analysis.")


#def save_analysis_result(self, out_name):
#    if out_name == None: out_name = result + '.pickle'
        
#    if result in list(self.analyses.keys()):
#        with open(out_name, 'wb') as handle:
#            pickle.dump(result, handle, protocol=2)

#        print(f"{analysis} analysis has been saved as {out_name}!")

    

