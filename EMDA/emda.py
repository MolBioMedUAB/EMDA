
from dataclasses import dataclass

# load MDAnalysis' universe class
from MDAnalysis import Universe

# load internal EMDA classes and functions
from selection import selection
import adders

# load custom exceptions
from exceptions import NotAvailableOptionError



class EMDA:

    def __init__(self, parameters, trajectory):

        self.parameters = parameters
        self.trajectory = trajectory
        self._universe  = Universe(parameters, trajectory)
        print('Trajectory has been loaded!')
        self.measures   = {}
        self.selections = {}

    @dataclass
    class Measure:
        name    : str
        type    : str
        sel     : list
        options : dict
        result  : list = []

    # function for creating selections (AtomGroups) as a dictionary inside EMDA class
    def select(self, name, sel_input, sel_type=None, no_backbone=False, return_atomic_sel_string=False):
        self.selections[name] = selection(self._universe, sel_input, sel_type=sel_type, no_backbone=no_backbone, return_atomic_sel_string=return_atomic_sel_string)


    # functions for creating Measures classes
    def add_distance(self, name, sel1, sel2, type='min'):
        """
        DESCRIPTION 
            Add a distance measure to be computed
        """
        #to_add = adders.distance(name, self.selections[sel1], self.selections[sel2], type=type)

        self.measures[name] = self.Measure(
            name, type, sel, options = adders.distance(name, self.selections[sel1], self.selections[sel2], type=type)
        )

        #self.measures[name] = self.Measure(
        #        name    = to_add['name'],
        #        type    = to_add['type'],
        #        sel     = to_add['sel'],
        #        options = to_add['options'],
        #        result  = []
        #    )
        

    def add_angle(self, name, sel1, sel2, sel3, units='deg', domain=360):
        to_add = adders.angle(name, self.selections[sel1], self.selections[sel2], self.selections[sel3], units=units, domain=domain)

        self.measures[name] = self.Measure(
                name    = to_add['name'],
                type    = to_add['type'],
                sel     = to_add['sel'],
                options = to_add['options'],
                result  = []
            )
        
    def add_dihedral(self, name, sel1, sel2, sel3, sel4, units='deg', domain=360):
        to_add = adders.angle(name, self.selections[sel1], self.selections[sel2], self.selections[sel3], self.selections[sel4], units=units, domain=domain)

        self.measures[name] = self.Measure(
                name    = to_add['name'],
                type    = to_add['type'],
                sel     = to_add['sel'],
                options = to_add['options'],
                result  = []
            )

    def add_planar_angle(self, name, sel1, sel2, units='deg', domain=360):
        """
        DESCRIPTION 
            Add a planar angle measure to be computed 
        """
        to_add = adders.planar_angle(name, self.selections[sel1], self.selections[sel2], units=units, domain=domain)

        self.measures[name] = self.Measure(
                name    = to_add['name'],
                type    = to_add['type'],
                sel     = to_add['sel'],
                options = to_add['options'],
                result  = []
            )


    def add_contacts(self, name, sel, sel_env, interactions="all", include_WAT=False, out_format='new', measure_distances=True):
        """
        DESCRIPTION
            Add contacts measure to be computed
        """

        to_add = adders.contacts(self._universe, name, sel, sel_env, interactions=interactions, include_WAT=include_WAT, out_format=out_format, measure_distances=measure_distances)

        self.measures[name] = self.Measure(
            name    = to_add['name'],
            type    = to_add['type'],
            sel     = to_add['sel'],
            options = to_add['options'],
            result  = []
        )


    def add_RMSD(self, name, sel, ref=None, superposition=True):
        to_add = adders.RMSD(self._universe, name, sel, ref=ref, superposition=superposition)

        self.measures[name] = self.Measure(
            name    = to_add['name'],
            type    = to_add['type'],
            sel     = to_add['sel'],
            options = to_add['options'],
            result  = []
        )

    def add_distWATbridge(self, name, sel1, sel2, sel1_rad=3, sel2_rad=3):
        to_add = adders.distWATbridge(self._universe, name, sel1, sel2, sel1_rad=sel1_rad, sel2_rad=sel2_rad)

        self.measures[name] = self.Measure(
            name    = to_add['name'],
            type    = to_add['type'],
            sel     = to_add['sel'],
            options = to_add['options'],
            result  = []
        )

    def add_pKa(self, name, excluded_ions=["Na+", "Cl-"], pka_ref='neutral', pdb_folder='.pka', keep_pdb=False, keep_pka=False):
        print('WARNING!!:')
        print('  Take into account that pKa calculation with PROpKA requires')
        print('  the generation of a PDB for each analysed frame.')
        print('  This can be very time consuming for long trajectory, so consider')
        print('  using the step option in run_measure() to reduce')
        print('  the amount of analysed structures.')



        to_add = adders.pKa(self._universe, name, excluded_ions=excluded_ions, pka_ref=pka_ref, pdb_folder=pdb_folder, keep_pdb=keep_pdb, keep_pka=keep_pka)

        self.measures[name] = self.Measure(
            name    = to_add['name'],
            type    = to_add['type'],
            sel     = to_add['sel'],
            options = to_add['options'],
            result  = []
        )
    
    ##### END ADDING AVAILABLE MEASURES #####

    def remove_measurement(self, name):
        if isinstance(name, str): name = [name] 

        for measurement_index in range(len(self.measurements)):
            for name_ in name:
                if self.measurements[measurement_index]['name'] == name_:
                    self.measurements.pop(measurement_index)




emda = EMDA('../example/parameters.prmtop', '../example/trajectory.nc')
emda.select('first_resids', [1, 2, 3], sel_type='res_num')
emda.select('second_resids', [4, 5, 6], sel_type='res_num')

emda.add_distance('dist_first_second', 'first_resids', 'second_resids')

print(emda.measures)
