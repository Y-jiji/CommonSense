import os

def change_to_absolute(path, conf):
    if path[0] != "/":
        return os.path.join(conf["root path"], path)
    else:
        return path