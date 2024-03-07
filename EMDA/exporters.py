from MDAnalysis import Writer
from MDAnalysis.exceptions import *
from tqdm.autonotebook import tqdm

from .exceptions import *



def export_trajectory(self, variant, replica, out_name : str = None, split_in : int = 1, format : str = 'xtc', start : int = 0, end : int = -1, step : int = 1, selection = None):
    """
    DESCRIPTION:
        Function for exporting trajectories. Useful for exporting trajectories with transformations applied (such as wrap, unwrap or no jump).

    ARGUMENTS:
        - variant :     variant name to export
        - replica :     replica to export
        - out_name:     name for the output file. If * in the name, the number of split (if > 1) will be used. Conversely, the number of the split will be appended at the end.
        - split_in :    number of fragments to split the trajectory
        - format :      format of the output trajectory. It will be guessed from the output name. All files accepted  by MDAnalysis can be used (https://docs.mdanalysis.org/1.0.0/documentation_pages/coordinates/init.html#id2).
        - start :       first frame to be saved
        - end :         last frame to be saved
        - step :        number of frames to be stepped
    """

    # Check that variant and replica exist
    if variant not in list(self.universe.keys()):
        raise KeyError("Variant is not available")
    else :
        if replica not in list(self.universe[variant].keys()):
            raise KeyError(f"Replica is not available for variant {variant}")
        
    # load universe form EMDA
    universe = self.universe[variant][replica]

    # check/create out_name
    if out_name == None:
        if split_in > 1:
            out_name == f"md_{variant}_{replica}_*.{format}"

        elif split_in == 1:
            out_name == f"md_{variant}_{replica}.{format}"

    compatible_formats = ['dcd', 'xtc', 'trr', 'xyz', 'nc', 'pdb', 'crd', 'trz', 'mol2', 'coor', 'namdbin', 'in']

    if '.' in out_name:
        if out_name.split('.')[-1] not in compatible_formats:
            if split_in == 1:
                out_name = '.'.join([out_name, format])
            else :
                if '*' in out_name:
                    out_name = '.'.join([out_name, format])
                else :
                    out_name = '.'.join([out_name + '_*', format])
        
        else :
            if split_in >= 1:
                if '*' not in out_name:
                    out_name = '.'.join(out_name.split('.')[:-1]) + '_*' + out_name.split('.')[-1]

    else :
        if split_in >= 1:
            out_name = out_name + '_*' + format
        elif split_in == 1:
            out_name = out_name + format

    # Save trajectory 
    if split_in == 1:

        with Writer(out_name, universe.atoms.n_atoms) as W:
            for ts in tqdm(universe[start:end:step]):
                W.write(universe.atoms)

    else :
        split_lenght = len(universe[start:end])/split_in

        for split in tqdm(range(0,split_in), desc="Split"):
            if split == split_in-1:

                if selection == None:
                    with Writer(out_name, universe.atoms.n_atoms) as W:
                        print(f"Printing from frame {split_lenght*split} to last.")
                        for ts in tqdm(universe[split_lenght*split:-1:step]):
                            W.write(universe.atoms)

                else :
                    with Writer(out_name, universe.select_atoms(self.selections[selection]).atoms.n_atoms) as W:
                        print(f"Printing from frame {split_lenght*split} to last.")
                        for ts in tqdm(universe[split_lenght*split:-1:step]):
                            W.write(universe.select_atoms(self.selections[selection]).atoms)

            else :
                if selection == None:
                    with Writer(out_name, universe.atoms.n_atoms) as W:
                        print(f"Printing from frame {split_lenght*split} to {split_lenght*(split+1)}.")
                        for ts in tqdm(universe[split_lenght*split:split_lenght*(split+1):step]):
                            W.write(universe.atoms)

                else :
                    with Writer(out_name, universe.select_atoms(self.selections[selection]).atoms.n_atoms) as W:
                        print(f"Printing from frame {split_lenght*split} to {split_lenght*(split+1)}.")
                        for ts in tqdm(universe[split_lenght*split:split_lenght*(split+1):step]):
                            W.write(universe.select_atoms(self.selections[selection]).atoms)

            






