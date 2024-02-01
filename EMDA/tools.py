from os import path, mkdir, chdir

def check_folder(folder):

    folder = folder.split("/")
    current_path = path.abspath(".")

    for f in folder:
        if f != "":
            if not path.isdir(f):
                mkdir(f)

            chdir(f)

    chdir(current_path)

def get_most_frequent(list_):
    """
    Copied from https://www.geeksforgeeks.org/python-find-most-frequent-element-in-a-list/
    """
    counter = 0
    num = list_[0]
     
    for i in list_:
        curr_frequency = list_.count(i)
        if(curr_frequency> counter):
            counter = curr_frequency
            num = i
 
    return num


def read_analysis(self, analysis_filename, analysis=None):
    import pickle

    if analysis == None: analysis = '_'.join(analysis_filename.split('.')[:-1])

    with open(analysis_filename, 'rb') as handle:
        self.analyses[analysis] = pickle.load(handle)