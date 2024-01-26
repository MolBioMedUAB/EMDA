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
