""" common to all projects: utility tool for loading configurations """

"""In json files, all paths are given in relative path from the 'root path'.
   We do so for portability: When moving to another computer, you only need to 
   copy everything and change the root path in config. 
   In the loading, these paths are translated to absolute paths to prevent errors in programs."""

import json, os, random, copy, pathlib

def load_config(config_file):
    with open(config_file) as fin:
        config = json.load(fin)
    config_ori = copy.deepcopy(config)

    # translate relative paths to absolute paths
    for key, value in config.items():
        if "file" in key and value[0] != "/":
            config[key] = os.path.join(config["root path"], value)
            #mkdir if the request path does not exist
            parent_path = pathlib.Path(config[key]).parent
            parent_path.mkdir(parents=True, exist_ok=True)

    # generate random seeds
    if config.get("seed", "") == "RANDOM":
        config["seed"] = random.getrandbits(32)
        config_ori["seed"] = random.getrandbits(32)

    # write back random seeds
    with open(config_file, "w") as fout:
        json.dump(config_ori, fout, indent=4)

    return config
