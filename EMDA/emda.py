
from dataclasses import dataclass

# load MDAnalysis' universe class
from MDAnalysis import Universe
#from MDAnalysis.core.groups import AtomGroup

# load internal EMDA classes and functions
from .selection import selection
from .adders import *
from .runners import *

# load custom exceptions
from .exceptions import EmptyMeasuresError

from tqdm.autonotebook import tqdm

#from metaclass import add_adders
#class EMDA(metaclass=add_adders):

class EMDA:

    def __init__(self, parameters, trajectory):

        self.parameters = parameters
        self.trajectory = trajectory
        self.universe  = Universe(parameters, trajectory)
        print('Trajectory has been loaded!')
        self.measures   = {}
        self.selections = {}
        
        # Automatically add all functions from adders.py
        external_functions = [func for func in globals() if callable(globals()[func]) and func.startswith("add_")]
        for func_name in external_functions:
            setattr(EMDA, func_name, globals()[func_name])

        #self.__add_adders()

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
            print_ += f"\tName: {self.name}\n"
            print_ += f"\tType: {self.type}\n"
            print_ += f"\tSel: {self.sel}\n"
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


    def run(self, exclude=None, step=1, start=1, end=-1,):
        """ 
        DESCRIPTION:
            Run all the measurements configured in self.measures

        OPTIONS:
            - exclude:  skip measures with the given name 
            - step:     Frames to jump during the analysis. Default is 1, so all the trajectory will be analysed.
            - start:    First frame to start the analysis. Default is 0.
            - end:      Last frame to analyse (included). Default is last frame of trajectory
        """
        
        if len(self.measures) == 0:
            raise EmptyMeasuresError
        
        if end == -1: end = len(self.universe.trajectory) # so if no last frame given, last in trajectory is chosen.

        if exclude == None: exclude = []

        # trajectory cycle
        for ts in tqdm(self.universe.trajectory[start-1:end:step], desc='Measuring', unit='Frame'):
            # measures cycle

            for measure in (set(self.measures.keys()) - set(exclude)):
                """
                TO-DO:
                    Try if there is an automatic way to run the proper runner depending on the type set in the Measure. \
                    All runners should be named like run_$type_of_measure, so maybe there is a chance to automate it.

                NOTES:
                    'continue' allows the code to not test all if statements once the executed one is done. Maybe since \
                    if/elif are used, there is no need to explicitly set continue.
                """

                if self.measures[measure].type == 'distance':
                    run_distance(self.measures[measure])
                    continue

                elif self.measures[measure].type == 'angle':
                    run_angle(self.measures[measure])
                    continue

                elif self.measures[measure].type == 'dihedral':
                    run_dihedral(self.measures[measure])
                    continue

                elif self.measures[measure].type == 'planar_angle':
                    run_planar_angle(self.measures[measure])
                    continue

                elif self.measures[measure].type == 'contacts':
                    run_contacts(self.measures[measure])
                    continue

                elif self.measures[measure].type == 'RMSD':
                    run_RMSD(self.measures[measure])
                    continue
                
                elif self.measures[measure].type == 'distWATbridge':
                    run_distWATbridge(self.measures[measure])
                    continue

                elif self.measures[measure].type == 'pka':
                    if len(self.measures[measure].result) == 0:
                        check_folder(self.measures[measure].options['pdb_folder'])

                    run_pka(self.measures[measure])
                    continue

                else :
                    print(f"The {measure} type is not available. If you need, you can create it.")

    def analyse():
        pass

    
    

