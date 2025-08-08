import json, glob, os
from typing import Callable

"""
    example key dictionary
    keys = {                              custom sum                      
        "points": ["avg", "max", "min", ("custom", 0, lambda x, y, dik: x+y), 
        ("avg", "key", lambda dik: dik["points"][0], "if", lambda dik: dik["time"] < 1000),   average if
        "time": "sum",
        "valid": "sum",
    }
"""

init_values = {"sum": 0, "count": 0, "min": float('inf'), "max": float('-inf'), "avg": (0, 0), "set": []}
operations = {"sum": lambda x, y: x + y, 
              "count": lambda x, y: x + 1,
    "avg": lambda x, y: (x[0] + y, x[1] + 1), "min": min, "max": max,
    "set": lambda x, y: x + [y]}
post_ops = {"avg": lambda x: x[0] / x[1] if x[1] > 0 else None}

def initialize_item(item):
    if isinstance(item, str):
        return init_values[item]
    elif isinstance(item, tuple):
        if item[0] == "custom":
            return item[1]
        else:
            return init_values[item[0]]
    else:
        print("Error: Each item in key dictionary must be a string or a tuple.")
        print(item)
        raise ValueError
    
def update_value(key, mode, dik, result):
    if isinstance(mode, str):
        return operations[mode](result, dik[key])
    elif isinstance(mode, tuple):
        if mode[0] == "custom":
            return mode[2](result, dik[key], dik)
        else:
            # if condition not satisfied
            if "if" in mode and not mode[mode.index("if") + 1](dik):
                return result
            if "key" in mode:
                return operations[mode[0]](result, mode[mode.index("key") + 1](dik))
            return operations[mode[0]](result, dik[key])

def post_process(mode, result):
    if isinstance(mode, str):
        return post_ops.get(mode, lambda x: x)(result)
    elif isinstance(mode, tuple):
        if mode[0] == "custom" and len(mode) > 3:
            return mode[3](result)
        else:
            return post_ops.get(mode[0], lambda x: x)(result)

def aggregate_result(file_pattern, keys : dict, remove_error = False):
    result = keys.copy()
    for key, mode in result.items():
        if isinstance(mode, str) or isinstance(mode, tuple):
            result[key] = initialize_item(mode)
        elif isinstance(mode, list):
            result[key] = []
            for v in mode:
                result[key].append(initialize_item(v))
        else:
            print("Error: Each item in key dictionary must be a string or a list.")
            return {}

    for file in glob.glob(file_pattern):
        try:
            with open(file, 'r') as f:
                f_dict = json.load(f)
        except Exception as e:
            print(f"File: {file}\tError: {e}")
            if remove_error:
                os.remove(file)
            continue

        for key, mode in keys.items():
            if isinstance(mode, str) or isinstance(mode, tuple):
                result[key] = update_value(key, mode, f_dict, result[key])
            elif isinstance(mode, list):
                for i, v in enumerate(mode):  
                    try:
                        result[key][i] = update_value(key, v, f_dict, result[key][i])
                    except Exception as e:
                        print(f"File: {file}\tKey: {key}\tMode: {mode}\nError: {e}")
                        raise e

    for key, mode in keys.items():
        if isinstance(mode, str) or isinstance(mode, tuple):
            result[key] = post_process(mode, result[key])
        elif isinstance(mode, list):
            for i, v in enumerate(mode):
                result[key][i] = post_process(v, result[key][i])
                     
    return result