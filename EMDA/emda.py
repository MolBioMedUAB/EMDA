from dataclasses import dataclass, field, is_dataclass
import pickle

from typing import Literal

# from typing import Union
# from types import NoneType

# load MDAnalysis' universe class
from MDAnalysis import Universe

# from MDAnalysis.core.groups import AtomGroup

# load internal EMDA classes and functions
from .selection import selection, parse_selection
from .adders import *
from .runners import *
from .analysers import *
from .plotters import *

#from .tools import get_dictionary_structure

# load custom exceptions
from .exceptions import EmptyMeasuresError, NotAvailableVariantError

from tqdm.autonotebook import tqdm

# from metaclass import add_adders
# class EMDA(metaclass=add_adders):

"""
IDEAS:
    - [] Make parameters and trajectory not only to strings but to dict with mutants and replicas.
"""


class EMDA:

    def __init__(self, parameters, trajectory=None, variant_name=None):
        """
        DESCRIPTION:
            Function to initialise the EMDA class by loading the parameters and trajectory as a MDAnalysis universe and loading adders, analysers and plotters as internal methods.

        ATTRIBUTES:
            - parameters:   name of the parameters and topology file as a dict containing the name of variant(s) and replica(s) 
          X - trajectory:   name or list of names of the trajectory file(s)
            - universe:     MDAnalysis Universe object containing the parameters and trajectory set as input of the class
            - selections:   Dictionary containing as key the name (ID) of a selection and the MDAnalysis AtomGroup object as value
            - measures:     Dictionary containing as key the name (ID) of a measure and the EMDA's Measure object as value
            - analyses:     Dictionary containing as key the name (ID) of an analysis and the EMDA's Analysis object as value


        METHODS:
            - add_*:        Adders loaded from adders.py file. The available adders and their description and usage can be printed using the print_available_adders EMDA's method
            - analyse_*:    Analysers loaded from analysers.py file.
            - plot_*:       Analysers loaded from plotters.py file. External plotter (indicated by using ext_ as function's name prefix) are not loaded.
        """

        #self.parameters = parameters
        #self.trajectory = trajectory

        # Check structure of variant(s) and trajectory(ies)
        if variant_name == None:
            variant_name = "V1"

        if isinstance(parameters, str) and trajectory == None:
            self.universe = { variant_name : 
                             { "R1" : Universe(parameters, trajectory) }
                            }
            self.parameters = { variant_name : parameters } 
            self.__variants = 1
            self.__replicas = 1

        elif isinstance(parameters, str) and isinstance(trajectory, (str, list)):
            self.universe = { variant_name : 
                             { "R1" : Universe(parameters, trajectory) }
                            }
            self.parameters = { variant_name : parameters } 
            self.__variants = 1
            self.__replicas = 1
        
        print("Trajectory has been loaded!")

        self.selections = { 'all' : [] }
        self.measures = {}
        self.analyses = {}

        # Automatically add all imported functions from adders.py and from analysers.py as EMDA methods
        external_functions = [
            func
            for func in globals()
            if callable(globals()[func]) and func.startswith("add_")
        ]
        external_functions += [
            func
            for func in globals()
            if callable(globals()[func]) and func.startswith("analyse_")
        ]
        external_functions += [
            func
            for func in globals()
            if callable(globals()[func]) and func.startswith("plot_")
        ]

        for func_name in external_functions:
            setattr(EMDA, func_name, globals()[func_name])


    @dataclass
    class Measure:
        """
        DESCRIPTION:
            Dataclass that stores calculated measures and related attributes.

        ATTRIBUTES:
            - name:     Name (ID) of the measure
            - type:     Type of the measure (distance, angle, dihedral, planar_angle, RMSD, and contacts are currently available)
            - sel:      Selections related to the measure as AtomGroups
            - options:  Empty dictionary containing different options to set the measure calculation
            - result:   List containing the measured results.

        METHODS:
            - plot:     Creates a simple plot of the calculated measures. Only available for distance, angle, dihedral, planar_angle, and RMSD types
        """

        name: str
        type: str
        sel: list
        options: dict
        result: list

        def __str__(self) -> str:
            if len(self.result) == 0:
                status = "Not calculated"
            else:
                status = "Calculated"

            print_ = f"Measure dataclass with:\n"
            print_ += f"\tName:   {self.name}\n"
            print_ += f"\tType:   {self.type}\n"
            print_ += f"\tSel:    {self.sel}\n"
            print_ += f"\tStatus: {status}\n"

            return print_

        def __repr__(self):
            return self.__str__()

        def plot(self):
            """
            DESCRIPTION:
                Measure's method to plot the stored values in the result attribute for distance, angle, dihedral, RMSD and planar_angle types.
            """

            if self.type in ("distance", "angle", "dihedral", "RMSD", "planar_angle"):
                units = {
                    "distance": "(Å)",
                    "angle": "(º)",
                    "dihedral": "(º)",
                    "planar_angle": "(º)",
                    "RMSD": "(Å)",
                }

                plt.plot(self.result)
                plt.ylabel(
                    " ".join(self.type.split("_")).capitalize() + " " + units[self.type]
                )
                plt.xlabel("Frame")

                plt.show()
                plt.close()
            else:
                print(f"This method is still not available for {type} measures.")

    @dataclass
    class Analysis:
        """
        DESCRIPTION:
            Dataclass that stores calculated analyses and related attributes.

        ATTRIBUTES:
            - name:             Name (ID) of the measure
            - type:             Type of the measure (distance, angle, dihedral, planar_angle, RMSD, and contacts are currently available)
            - measure_name:     Selections related to the measure as AtomGroups
            - options:          Empty dictionary containing different options to set the measure calculation
            - result:           List containing the measured results.

        METHODS:
            - plot:     Creates a simple plot of the calculated measures. Only available for distance, angle, dihedral, planar_angle, and RMSD types
        """

        name: str
        type: str
        measure_name: str
        result: list
        options: dict = field(default_factory={})
        # options : dict
        # mode : Union[str, NoneType] = None

        def __str__(self) -> str:
            if len(self.result) == 0:
                status = "Not calculated"
            else:
                status = "Calculated"

            print_ = f"Analysis dataclass with:\n"
            print_ += f"\tName:   {self.name}\n"
            print_ += f"\tType:   {self.type}\n"
            print_ += f"\tRelated mesure:    {self.measure_name}\n"
            print_ += f"\tStatus: {status}\n"

            return print_

        def __repr__(self):
            return self.__str__()


    def load_variant(self, parameters, trajectory, name=None):
        """
        DESCRIPTION:
            Method that allows adding one new variant to the EMDA class. Its key in the EMDA's universe attr's dictionary is automatically given as "V" and the number of variant in __variants attr.
        """

        self.__variants += 1
        self.__replicas += 1

        if name == None:
            new_variant  = f"V{self.__variants}"
        else :
            new_variant = name

        self.universe[new_variant]   = {"R1" : Universe(parameters, trajectory)}
        self.parameters[new_variant] = parameters

        print(f"{new_variant} variant has been loaded!")


    def load_trajectory(self, trajectory, variant_name='last'):
        """
        DESCRIPTION:
            Method that allows adding one more replica to a pre-existing variant in th EMDA class.
        """

        if variant_name == 'last':
            variant_name = list(self.universe.keys())[-1]

        else :
            if variant_name not in list(self.parameters.keys()):
                raise NotAvailableVariantError


        new_replica = int(max(list(self.universe[variant_name].keys()))[1:]) + 1

        self.universe[variant_name][f"R{new_replica}"] = Universe(self.parameters[variant_name], trajectory)

        
        print(f"A new replica has been loaded to variant {variant_name}!")


    # function for creating selections (AtomGroups) as a dictionary inside EMDA class
    def select(
        self,
        name,
        sel_input=None,
        sel_type=None,
        no_backbone=False,
        #variant=None
    ):
        """
        DESCRIPTION:
            Method for creating selections inside the selections attribute of EMDA's class.
        """

        if sel_input == None:
            sel_input = name

        # Creates the selection for all variants
        
        self.selections[name] = parse_selection(sel_input=sel_input, sel_type=sel_type, no_backbone=no_backbone)

        
        # Creates selection for each variant
        #if variant == None:
        #    for variant_ in list(self.universe.keys()):
        #        try :
        #            self.selections[variant_][name] = parse_selection(sel_input=sel_input, sel_type=sel_type, no_backbone=no_backbone)
        #        except KeyError:
        #            self.selections[variant_] = {}
        #            self.selections[variant_][name] = parse_selection(sel_input=sel_input, sel_type=sel_type, no_backbone=no_backbone)
        #
        #else :
        #    if variant in list(self.universe.keys()):
        #        self.selections[variant][name] = parse_selection(sel_input=sel_input, sel_type=sel_type, no_backbone=no_backbone)
        #    else :
        #        raise NotAvailableVariantError(variant)
        #    
        #self.selections['all'].append(name)


    def print_available_adders(self):
        """
        DESCRIPTION:
            Prints the DESCRIPTION and USAGE sections of the available adders
        """

        import inspect

        adder_methods = [method for method in dir(self) if method.startswith("add_")]

        for method_name in adder_methods:
            method = getattr(self, method_name)

            docstring = inspect.getdoc(method).split("\n")

            # Extract DESCRIPTION and USAGE sections from docstring
            description_section = ""
            usage_section = ""
            for l in range(len(docstring)):
                if "DESCRIPTION" in docstring[l]:
                    description_section += docstring[l] + "\n"
                    l_ = l
                    while True:
                        l_ += 1
                        if docstring[l_] == "":
                            break
                        else:
                            description_section += docstring[l_] + "\n"

                if "USAGE" in docstring[l]:
                    usage_section += docstring[l] + "\n"
                    l_ = l
                    while True:
                        l_ += 1
                        if docstring[l_] == "":
                            break
                        else:
                            usage_section += docstring[l_] + "\n"

            if description_section == "":
                description_section = "No description available.\n"
            if usage_section == "":
                usage_section = "No usage available.\n"

            print(
                f"\033[1mMethod: {method_name}\033[0m\nHelp:\n{description_section}\n{usage_section}\n"
            )

        print(
            "Use '>>> help(EMDA.add_***)' to get the complete information of an adder."
        )

    def run(
        self,
        exclude=None,
        run_only=None,
        recalculate=False,
        step=1,
        start=1,
        end=-1,
    ): 
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
        if end == -1:
            end = len(self.universe.trajectory)

        # Convert exclude to list to append precalculated measures if recalculate is False.
        ## If it is True, set measure's result as empty list, so it is overwritten.
        ## If a measure or list of measures is given, their result list will be reset as [].
        if exclude == None:
            exclude = []
        elif isinstance(exclude, str):
            exclude = [exclude]
        ###
            
        # Checks if there is contents in results and resets the list 
        if isinstance(recalculate, str):
            recalculate = [recalculate]

        for measure in list(self.measures.keys()):
            if len(self.measures[measure].result) > 0:
                if isinstance(recalculate, bool):
                    if recalculate:
                        self.measures[measure].result = []
                    elif not recalculate:
                        exclude.append(measure)

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
        elif run_only == None:
            measures = set(set(self.measures.keys()) - set(exclude))

        # trajectory cycle
        first_cycle = True
        for ts in tqdm(
            self.universe.trajectory[start - 1 : end : step],
            desc="Measuring",
            unit="Frame",
        ):

            # measures cycle
            for measure in measures:
                """
                TO-DO:
                    Look for an automatic way to run the proper runner depending on the type set in the Measure. \
                    All runners should be named like run_${type_of_measure}, so maybe there is a chance to automate it.
                """

                if self.measures[measure].type == "distance":
                    run_distance(self.measures[measure])

                elif self.measures[measure].type == "angle":
                    run_angle(self.measures[measure])

                elif self.measures[measure].type == "dihedral":
                    run_dihedral(self.measures[measure])

                elif self.measures[measure].type == "planar_angle":
                    run_planar_angle(self.measures[measure])

                elif self.measures[measure].type == "contacts":
                    run_contacts(self.measures[measure])

                elif self.measures[measure].type == "RMSD":
                    run_RMSD(self.measures[measure])

                elif self.measures[measure].type == "distWATbridge":
                    run_distWATbridge(self.measures[measure])

                elif self.measures[measure].type == "pka":
                    if first_cycle:
                        check_folder(self.measures[measure].options["pdb_folder"])

                    run_pka(self.measures[measure])

                else:
                    if first_cycle:
                        print(
                            f"The {measure} type is not available. If you need, you can create it."
                        )

            if first_cycle:
                first_cycle = False

    #def save_result(self, name, out_name=None):
    #    """
    #    DESCRIPTION:
    #        EMDA's method for saving into a pickle file the result attribute of a Measure class or an Analysis class.
    #    """
    #
    #    if out_name == None:
    #        out_name = name + ".pickle"
    #
    #    if name in list(self.analyses.keys()):
    #        with open(
    #            ".".join(out_name.split(".")[:-1])
    #            + "_analysis."
    #            + out_name.split(".")[-1],
    #            "wb",
    #        ) as handle:
    #            pickle.dump(self.analyses[name].result, handle, protocol=2)
    #
    #        print(f"{name} analysis' result has been saved as {out_name}!")
    #
    #    elif name in list(self.measures.keys()):
    #        with open(
    #            ".".join(out_name.split(".")[:-1])
    #            + "_measure."
    #            + out_name.split(".")[-1],
    #            "wb",
    #        ) as handle:
    #            pickle.dump(self.measures[name].result, handle, protocol=2)
    #
    #        print(f"{name} analysis' result has been saved as {out_name}!")
    #
    #    else:
    #        raise KeyError(f"{name} is not an available measure nor analysis.")
    #
    #def read_result(
    #    self,
    #    filename,
    #    name,
    #    type: Literal["d", "dataclass", "a", "analysis", "m", "measure"] = "dataclass",
    #):
    #    """
    #    DESCRIPTION:
    #        EMDA's method for reading a pickle file containing the result attribute of a precreated Measure class or an Analysis class.
    #    """
    #
    #    if is_dataclass(name) and type.lower() in ["d", "dataclass"]:
    #        with open(filename, "rb") as handle:
    #            name.result = pickle.load(handle)
    #
    #    else:
    #        try:
    #            if type.lower() in ["m", "measure"]:
    #                with open(filename, "rb") as handle:
    #                    self.measures[name].result = pickle.load(handle)
    #
    #            if type.lower() in ["a", "analysis"]:
    #                with open(filename, "rb") as handle:
    #                    self.analyses[name].result = pickle.load(handle)
    #
    #        except KeyError:
    #            raise KeyError(f"{name} is not an available measure nor analysis.")
