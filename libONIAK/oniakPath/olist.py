import os, re

""" List all files under the pattern. """
def list_pattern(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for file in files:
            full_path = os.path.join(root, file)
            if re.search(pattern, full_path):
                yield full_path
