from oniakIO import odats
from glob import glob
import os

# check how many files follow the given pattern and has the given size
# if delete is true, files of wrong sizes will be deleted
def check_files(pattern, size, delete=False):
    cnt = 0
    for file in glob(pattern):
        file_reader = odats.ONIAKReader(file)
        if file_reader.get_size() == size:
            cnt += 1
        elif delete:
           os.remove(file)
    return cnt