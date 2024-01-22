
from dataclasses import dataclass

# load MDAnalysis' universe class
from MDAnalysis import Universe
#from MDAnalysis.core.groups import AtomGroup

# load internal EMDA classes and functions
from selection import selection
from adders import *
# adders

# load custom exceptions
from exceptions import EmptyMeasuresError

from tqdm.autonotebook import tqdm

#from metaclass import add_adders
#class EMDA(metaclass=add_adders):

class EMDA:

    def __init__(self, parameters, trajectory):

        self.parameters = parameters
        self.trajectory = trajectory
        self.__universe  = Universe(parameters, trajectory)
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

    # function for creating selections (AtomGroups) as a dictionary inside EMDA class
    def select(self, name, sel_input, sel_type=None, no_backbone=False, return_atomic_sel_string=False):
        self.selections[name] = selection(self._universe, sel_input, sel_type=sel_type, no_backbone=no_backbone, return_atomic_sel_string=return_atomic_sel_string)


    def print_available_adders(self, all_description=False):

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
            Run all the measurements configured in self.measures. 

        OPTIONS:
            - exclude:  skip measures with the given name 
            - step:     Frames to jump during the analysis. Default is 1, so all the trajectory will be analysed.
            - start:    First frame to start the analysis. Default is 0.
            - end:      Last frame to analyse (included). Default is last frame of trajectory
        """
        
        if len(self.measures) == 0:
            raise EmptyMeasuresError
        
        if end == -1: end = len(self.__universe.trajectory) # so if no last frame given, last in trajectory is chosen.

        # trajectory cycle
        for ts in tqdm(self.__universe.trajectory[start, end, step], desc='Measuring', unit='Frame'):
            # measures cycle
            for measure in self.measures:
                pass







        pass

    def analyse():
        pass

    
    

emda = EMDA('../example/parameters.prmtop', '../example/trajectory.nc')
#emda.select('first_resids', [1, 2, 3], sel_type='res_num')
#emda.select('second_resids', [4, 5, 6], sel_type='res_num')

#emda.add_distance(name='dist_first_second', sel1='first_resids', sel2='second_resids', type='min')

print(emda.print_available_adders())

