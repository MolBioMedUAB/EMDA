from dataclasses import dataclass, field, is_dataclass
import pickle

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

from time import sleep

"""
IDEAS:
    - [] Make parameters and trajectory not only to strings but to dict with mutants and replicas.
"""


class EMDA:

    def __init__(self, parameters, trajectory=None, variant_name=None, load_in_memory : bool = False):
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
                             { "R1" : Universe(parameters, trajectory, in_memory=load_in_memory) }
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

        self.selections = {}
        self.measures = {}
        self.analyses = {}

        self.load_in_memory = load_in_memory

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
            print_ = f"Measure dataclass with:\n"
            print_ += f"\tName:   {self.name}\n"
            print_ += f"\tType:   {self.type}\n"
            print_ += f"\tSel:    {self.sel}\n"
            print_ += f"\tStatus: \n"

            #status = get_dictionary_structure(self.result, False)
            for variant in list(self.result.keys()):
                for replica in list(self.result[variant]):
                    if len(self.result[variant][replica]) > 0:
                        print_ += f"\t\t{variant}, {replica}: Calculated\n"
                        #status[variant][replica] = "Calculated"
                    else :
                        print_ += f"\t\t{variant}, {replica}: Not calculated\n"
                        #status[variant][replica] = "Not calculated"
 
            return print_

        def __repr__(self):
            return self.__str__()
        

        def plot(self, same_y : bool = True, same_x : bool = True, axis_label_everywhere=False, combine_replicas=False, out_name=None):
            plot_measure(self, measure_name=None, same_y=same_y, same_x=same_x, axis_label_everywhere=axis_label_everywhere, combine_replicas=combine_replicas, out_name=out_name)

    @dataclass
    class Analysis:
        """
        DESCRIPTION:
            Dataclass that stores calculated analyses and related attributes.

        ATTRIBUTES:
            - name:             Name (ID) of the analysis
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
            print_ = f"Analysis dataclass with:\n"
            print_ += f"\tName:   {self.name}\n"
            print_ += f"\tType:   {self.type}\n"
            print_ += f"\tRelated mesure:    {self.measure_name}\n"

            return print_

        def __repr__(self):
            return self.__str__()
        
        def plot(self, analysis_name=None, merge_replicas=False, percentage=False, error_bar=True, bar_width=None, width=None, errorbar_width=5 , width_per_replica=4, height_per_variant=4, title=None, same_y=True, same_x=True, axis_label_everywhere=False, residue_label_rotation=45, out_name=False):
            if self.type in ('value', 'NACs'):
                if bar_width == None:
                    bar_width = 0.1
                plot_NACs(self, analysis_name=analysis_name, merge_replicas=merge_replicas, percentage=percentage, error_bar=error_bar, bar_width=bar_width, width=width, title=title, out_name=out_name)
            
            elif self.type in ('contacts_amount') and self.options['mode'] in ('contacts'):
                plot_measure(self, measure_name=None, same_y=same_y, same_x=same_x, axis_label_everywhere=axis_label_everywhere, combine_replicas=merge_replicas, out_name=out_name)

            elif self.type in ("contacts_frequency") and self.options['mode'] in ('contacts'):
                if bar_width == None:
                    bar_width = 0.8
                plot_contacts_frequency(self, analysis_name=analysis_name, same_y=same_y, same_x=same_x, axis_label_everywhere=axis_label_everywhere, merge_replicas=merge_replicas, error_bar=error_bar, bar_width=bar_width, errorbar_width=errorbar_width, width_per_replica=width_per_replica, height_per_variant=height_per_variant, residue_label_rotation=residue_label_rotation, out_name=out_name)

            else :
                raise NotCompatibleAnalysisForPlotterError



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

        self.universe[new_variant]   = {"R1" : Universe(parameters, trajectory, in_memory=self.load_in_memory)}
        self.parameters[new_variant] = parameters

        # Adds new variant and replica to existing measures
        for measure in list(self.measures.keys()):
                self.measures[measure].result[new_variant] = {"R1" : []}

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

        self.universe[variant_name][f"R{new_replica}"] = Universe(self.parameters[variant_name], trajectory, in_memory=self.load_in_memory)

        # Adds new variant and replica to existing measures
        for measure in list(self.measures.keys()):
            self.measures[measure].result[variant_name][f"R{new_replica}"] = []
        
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
        verbose=False,
        sleep_time=0,
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
        
        else :
            pass

        # Creates a dict with the same structure as universe. If end is -1, the lenght of each trajectory is read and
        ## saved. If not -1, the length is checked and if traj is longer, the given end is saved.
        ends = get_dictionary_structure(self.universe, 0)
        for variant in list(self.universe.keys()):
            for replica, u in self.universe[variant].items():
                if end == -1 or end > len(u.trajectory):
                    ends[variant][replica] = len(u.trajectory)
                else :  
                    ends[variant][replica] = end

        # Crates a dict with the same structure as universe for starts and step.
        starts = get_dictionary_structure(self.universe, start-1)
        steps  = get_dictionary_structure(self.universe, step)

        # Convert exclude to list to append precalculated measures if recalculate is False.
        ## If it is True, set measure's result as empty list, so it is overwritten.
        ## If a measure or list of measures is given, their result list will be reset as [].
        if exclude == None:
            excludes = get_dictionary_structure(self.universe, [])
        elif isinstance(exclude, str):
            excludes = get_dictionary_structure(self.universe, [exclude])
        ###
            
        # Checks if there is contents in results and resets the list 
        recalculates = get_dictionary_structure(self.universe, [])
        if isinstance(recalculate, str):
            recalculates = get_dictionary_structure(self.universe, [recalculate])

        for measure in list(self.measures.keys()):
            for variant in list(recalculates.keys()):
                for replica in list(recalculates[variant].keys()):
                    if len(self.measures[measure].result[variant][replica]) > 0:
                        if isinstance(recalculate, bool):
                            if recalculate:
                                self.measures[measure].result[variant][replica] = []
                            elif not recalculate:
                                excludes[variant][replica].append(measure)

                        elif isinstance(recalculate, list):
                            if measure in recalculate:
                                self.measures[measure].result[variant][replica] = []

        # Check run_only. If run_only is used, the measure set will be the run_only list. Conversely, measures will be set as
        ## all measures except exclude (if none, it is converted to empty list in previous codeblock)
        if run_only != None:
            if isinstance(run_only, str):
                measures = get_dictionary_structure(self.universe, set([run_only]))
            elif isinstance(run_only, list):
                measures = get_dictionary_structure(self.universe, set([run_only]))
        elif run_only == None:
            measures = get_dictionary_structure(self.universe, set())
            for variant in list(self.universe.keys()):
                for replica in list(self.universe[variant].keys()):
                    measures[variant][replica] = set(set(self.measures.keys()) - set(excludes[variant][replica]))


        def run_measures(self, measures, variant, replica):

            #u = self.universe[variant][replica]

            # trajectory cycle
            first_cycle = True
            for ts in tqdm(self.universe[variant][replica].trajectory[starts[variant][replica] : ends[variant][replica] : steps[variant][replica]],
                           desc=f"Measuring variant {variant}, replica {replica}",
                           unit="Frame"
                        ):
                
                # measures cycle
                for measure in measures[variant][replica]:
                    """
                    TO-DO:
                        Look for an automatic way to run the proper runner depending on the type set in the Measure. \
                        All runners should be named like run_${type_of_measure}, so maybe there is a chance to automate it.
                    """

                    if self.measures[measure].type == "distance":
                        run_distance(self, self.measures[measure], variant=variant, replica=replica)

                    elif self.measures[measure].type == "angle":
                        run_angle(self, self.measures[measure], variant=variant, replica=replica)

                    elif self.measures[measure].type == "dihedral":
                        run_dihedral(self, self.measures[measure], variant=variant, replica=replica)

                    elif self.measures[measure].type == "planar_angle":
                        run_planar_angle(self, self.measures[measure], variant=variant, replica=replica)

                    elif self.measures[measure].type == "contacts":
                        run_contacts(self, self.measures[measure], variant=variant, replica=replica)

                    elif self.measures[measure].type == "protein_contacts":
                        run_protein_contacts(self, self.measures[measure], variant=variant, replica=replica)

                    elif self.measures[measure].type == "RMSD":
                        run_RMSD(self, self.measures[measure], variant=variant, replica=replica)

                    elif self.measures[measure].type == "distWATbridge":
                        run_distWATbridge(self, self.measures[measure], variant=variant, replica=replica)

                    elif self.measures[measure].type == "pka":
                        if first_cycle:
                            check_folder(self.measures[measure].options["pdb_folder"])

                        run_pka(self, self.measures[measure], variant=variant, replica=replica)

                    else:
                        if first_cycle:
                            print(
                                f"The {measure} type is not available. If you need, you can create it."
                            )

                if first_cycle:
                    first_cycle = False

                sleep(sleep_time)

            del ts, #u



        if len(self.universe) == 1:
            variant = list(self.universe.keys())[0]
            if len(self.universe[variant]) == 1:
                replica = list(self.universe[variant].keys())[0]
                run_measures(self, measures=measures, variant=variant, replica=replica)
            else :
                print('single variant, multireplica')
                for replica in tqdm(list(self.universe[variant].keys()),
                                    desc="Replica",
                                    unit="rep"
                                ):
                    run_measures(self, measures=measures, variant=variant, replica=replica)

        else :
            # variants cycle
            for variant in tqdm(list(self.universe.keys()), desc='Variants', unit='var'):
                if verbose:
                    print(f"Starting variant {variant} ")
                # replicas cycle
                r_num = 0
                for replica in list(self.universe[variant].keys()):
                    r_num += 1
                    if verbose:
                        print(f"Starting replica {replica} ({r_num} of {len(self.universe[variant].keys())})")
                    run_measures(self, measures=measures, variant=variant, replica=replica)




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
