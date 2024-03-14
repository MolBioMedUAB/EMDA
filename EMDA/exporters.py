from MDAnalysis import Writer
from MDAnalysis.exceptions import *
from tqdm.autonotebook import tqdm

from .exceptions import *

from .tools import check_folder


def export_frames_by_analysis(self, variant, replica, analysis_name, out_name=None, format='pdb', folder=None, selection = 'all'):
    """
    DESCRIPTION:
        Function for exporting trajectories. Useful for exporting trajectories with transformations applied (such as wrap, unwrap or no jump).

    ARGUMENTS:
        - variant :     variant name to export
        - replica :     replica to export
        - out_name:     name for the output file. Use * in the out_name to specify the position of the frame. If not present, it will be located before the extension.
        - format:       only pdb is available
        - folder:       subfolder where exported frames are saved

    TODO:
        - [ ] Add more formats
    """
    
    #compatible_formats = ['dcd', 'xtc', 'trr', 'xyz', 'nc', 'pdb', 'crd', 'trz', 'mol2', 'coor', 'namdbin', 'in']
    compatible_formats = ['pdb']


    # Check that variant and replica exist
    if variant not in list(self.universe.keys()):
        raise KeyError("Variant is not available")
    else :
        if replica not in list(self.universe[variant].keys()):
            raise KeyError(f"Replica is not available for variant {variant}")
        
    # load universe from EMDA
    universe = self.universe[variant][replica]

    if out_name == None:
        out_name = f"md_{variant}_{replica}_frame_*.pdb"

    # check out_name
    if '.' in out_name:
        if out_name.split('.')[-1] not in compatible_formats:
            if '*' in out_name:
                out_name = '.'.join([out_name, format])
            else :
                out_name = '.'.join([out_name + '_*', format])
        
        else :
            if '*' not in out_name:
                out_name = '.'.join(out_name.split('.')[:-1]) + '_*' + out_name.split('.')[-1]

    else :
        out_name = out_name + '_*' + format

    if folder != None:
        out_name = '/'.join([folder, out_name])
        check_folder(folder)


    if selection in list(self.selections):
        selection = self.selections[selection]

    for frame, result in enumerate(self.analyses[analysis_name].result[variant][replica]):

        if result:
            universe.trajectory[frame]
            to_write = universe.select_atoms(selection)
            to_write.write(out_name.replace('*', str(frame+1)))
            



def export_trajectory(self, variant, replica, out_name : str = None, split_in : int = 1, format : str = 'xtc', start : int = 0, end : int = -1, step : int = 1, selection = None, folder = None):
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
        
    # load universe from EMDA
    universe = self.universe[variant][replica]

    compatible_formats = ['dcd', 'xtc', 'trr', 'xyz', 'nc', 'pdb', 'crd', 'trz', 'mol2', 'coor', 'namdbin', 'in']

    # check/create out_name
    if out_name == None:
        if split_in > 1:
            out_name = f"md_{variant}_{replica}_*.{format}"

        elif split_in == 1:
            out_name = f"md_{variant}_{replica}.{format}"

    else :
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
            

    if folder != None:
        out_name = '/'.join([folder, out_name])
        check_folder(folder)
        
    # 
    if selection in list(self.selections):
        selection = self.selections[selection]

    if end == -1:
        end = len(universe.trajectory) + 1

    # Save trajectory 
    if split_in == 1:

        if selection == None:
            with Writer(out_name, universe.atoms.n_atoms) as W:
                for ts in tqdm(universe.trajectory[start:end:step]):
                    W.write(universe.atoms)

        else :
            with Writer(out_name, universe.select_atoms(selection).atoms.n_atoms) as W:
                for ts in tqdm(universe.trajectory[start:end:step], desc='Saving trajectory', unit='frame'):
                    W.write(universe.select_atoms(selection).atoms)


    else :
        split_length = int((len(universe.trajectory[start:end])+1)/split_in)

        for split in tqdm(range(0,split_in), desc="Split"):
            if split == split_in-1:

                if selection == None:
                    with Writer(out_name.replace('*', str(split+1)), universe.atoms.n_atoms) as W:
                        print(f"Printing from frame {split_length*split+1} to last ({end}).")
                        for ts in tqdm(universe.trajectory[split_length*split:end:step], desc='Saving trajectory', unit='frame'):
                            W.write(universe.atoms)

                else :
                    with Writer(out_name.replace('*', str(split+1)), universe.select_atoms(self.selections[selection]).atoms.n_atoms) as W:
                        print(f"Printing from frame {split_length*split+1} to last ({end}).")
                        for ts in tqdm(universe.trajectory[split_length*split:end:step], desc='Saving trajectory', unit='frame'):
                            W.write(universe.select_atoms(selection).atoms)

            else :
                if selection == None:
                    with Writer(out_name.replace('*', str(split+1)), universe.atoms.n_atoms) as W:
                        print(f"Printing from frame {split_length*split+1} to {split_length*(split+1)}.")
                        for ts in tqdm(universe.trajectory[split_length*split:split_length*(split+1):step], desc='Saving trajectory', unit='frame'):
                            W.write(universe.atoms)

                else :
                    with Writer(out_name.replace('*', str(split+1)), universe.select_atoms(self.selections[selection]).atoms.n_atoms) as W:
                        print(f"Printing from frame {split_length*split+1} to {split_length*(split+1)}.")
                        for ts in tqdm(universe.trajectory[split_length*split:split_length*(split+1):step], desc='Saving trajectory', unit='frame'):
                            W.write(universe.select_atoms(selection).atoms)

            






