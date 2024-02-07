from os import path, mkdir, chdir


def check_folder(folder):
    """
    DESCRIPTION:
        Function to check if a requested folder exists and move to it. If it does not exist, it is created and then the directory is changed.
    """

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
    DESCRIPTION:
        Function to obtain the most frequent value in a list and return it.

    SOURCE:
        Copied from https://www.geeksforgeeks.org/python-find-most-frequent-element-in-a-list/
    """
    counter = 0
    num = list_[0]

    for i in list_:
        curr_frequency = list_.count(i)
        if curr_frequency > counter:
            counter = curr_frequency
            num = i

    return num


# def read_analysis(self, analysis_filename, analysis=None):
#    import pickle

#    if analysis == None: analysis = '_'.join(analysis_filename.split('.')[:-1])

#    with open(analysis_filename, 'rb') as handle:
#        self.analyses[analysis] = pickle.load(handle)


def in_notebook() -> bool:
    """
    DESCRIPTION:
        Function for checking if the code is being executed within a Jupyter Notebook/IPython env. The returned value is of bool type.

    SOURCE
        https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    """

    try:
        from IPython import get_ipython

        if "IPKernelApp" not in get_ipython().config:  # pragma: no cover
            return False
    except ImportError:
        return False
    except AttributeError:
        return False
    return True
