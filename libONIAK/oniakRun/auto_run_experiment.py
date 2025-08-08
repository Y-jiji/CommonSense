import sys, json, os, subprocess, pathlib, copy
import argparse, time, random, psutil, glob, re
import numpy as np
import numbers
from string import Formatter

parser = argparse.ArgumentParser()
parser.add_argument("conf", nargs="+", help="conf jsons are required.")
parser.add_argument("--wait", action="store_true")
parser.add_argument("--waitload", default=0.3, type=float)
parser.add_argument("--waittime", default=90.0, type=float)
arguments = parser.parse_args()

if arguments.wait:
    # wait for the current jobs to finish
    while True:
        sysavg = os.getloadavg()[0]
        if sysavg < arguments.waitload:
            break
        else:
            time.sleep(arguments.waittime)


# This iterates and expands nested dictionaries in the config file
def expand_keys(dikt, key_list, stack):
    for key in dikt.keys():
        if isinstance(dikt[key], dict):
            expand_keys(dikt[key], key_list, stack + [key])
        else:
            key_list.append(stack + [key])

# read from nested dictionary. If key is a list, read until the end.
def read_dict(dikt, key):
    if isinstance(key, str):
        return dikt[key]
    a = dikt
    for k in key:
        a = a[k]
    return a


def write_dict(dikt, key, value):
    final_dikt = read_dict(dikt, key[:-1])
    final_dikt[key[-1]] = value


def run_conf(config_path):
    with open(config_path) as fin:
        config = json.load(fin)

    processes = set()
    max_processes = config["max process"]
    # max memory PERCENTAGE before a job is put on wait
    max_memory = config.get("max memory", 90)
    cnt = 0
    config_folder = pathlib.Path(config_path).parent
    config_name = pathlib.Path(config_path).name.split(".")[0]

    # conf is original config dictionary
    # conf_exp is config dictionary after expansion
    # key_list is a list of keys (each key is a list) in conf
    # level is current pointer at key_list
    # fill_dict is a flattened dictionary for string replacement
    def expand_recur(conf, conf_exp, key_list, level, fill_dict):    
        class LevelFormatter(Formatter):
            def get_value(self, key, args, kwds):
                if isinstance(key, str):
                    if level >= len(key_list) - 1 or key in [x[-1] for x in key_list[:level+1]]:
                        return kwds[key]
                    else:
                        return "*"
                else:
                    return Formatter.get_value(key, args, kwds)

        def replace_str(input):
            fmt = LevelFormatter()
            result = input.replace("_AUTO", "_{}".format(cnt))
            if "{" in input and "}" in input:
                try:
                    result = fmt.format(input, **fill_dict)
                except IndexError:
                    print(input)
                    raise IndexError("Error: Index Error in format string")
            if "ROOT" in result:
                root_path = pathlib.Path(conf["root path"])
                result = str(root_path / result[5:])  # must start with ROOT/
            return result

        nonlocal cnt
        while level < len(key_list):
            cur_key = key_list[level]
            cur_element = read_dict(conf, cur_key)
            if (
                isinstance(cur_element, list)
                and isinstance(cur_element[0], str)
                and cur_element[0][:4] == "ITER"
            ):
                iter_list = []
                # iterate over all options in the list
                if cur_element[0] == "ITER_ALL":
                    iter_list = cur_element[1:]
                # iterate over a range given by its end (start, or step)
                elif cur_element[0] == "ITER_RANGE":
                    iter_list = range(*cur_element[1:])
                # iterate over float numbers using np.linspace
                elif cur_element[0] == "ITER_LINSPACE":
                    iter_list = np.linspace(*cur_element[1:])

                # IMPORTANT: pay attention to where you use ITER_RANDOM
                # for example, if you iterate over both queries and seeds,
                # but want all queries to have the same seed
                # then you HAVE TO place "seed" before "query"
                # otherwise the seeds are independently generated for each

                elif cur_element[0] == "ITER_RANDOM":
                    iter_list = [random.getrandbits(32) for i in range(cur_element[1])]

                # iterate synchronously with an ABOVE key

                elif cur_element[0] == "ITER_WITH":
                    with_key = cur_element[1]
                    with_dict = cur_element[2]
                    assert with_key in fill_dict.keys()
                    with_key_value = str(fill_dict[with_key])
                    if with_key_value not in with_dict.keys():
                        with_key_value = "DEFAULT"
                    iter_list = with_dict[with_key_value]
                    # shorthand: list can be ignored if length is 1
                    if not isinstance(iter_list, list):
                        iter_list = [iter_list]

                # iterate over all file names with pattern
                # optional: length of iteration list
                elif cur_element[0] == "ITER_FILE":
                    pattern = cur_element[1]
                    iter_list = glob.glob(pattern)
                    if len(cur_element) == 3:
                        max_len = cur_element[2]
                        iter_list = iter_list[:max_len]
                # extract keyword from other string
                # can only extract integer numbers
                elif cur_element[0] == "ITER_EXTRACT":
                    with_key = cur_element[1]
                    extract_key = cur_element[2]
                    assert with_key in fill_dict.keys()
                    from_string = fill_dict[with_key]
                    re_pattern = "_" + extract_key + r"(\d*?)(_|\.odat|\.json)"
                    iter_list = [int(re.search(re_pattern, from_string).group(1))]
                else:
                    print("Error: ITER keyword unrecognized. Job NOT executed.")
                    return

                if len(iter_list) == 0:
                    print("Error: Empty Iter List. No job generated from this config")
                    return

                for val in iter_list:
                    write_dict(conf_exp, cur_key, val)
                    fill_dict[cur_key[-1]] = val
                    expand_recur(conf, conf_exp, key_list, level + 1, fill_dict)
                return
            elif (
                isinstance(cur_element, list)
                and isinstance(cur_element[0], str)
                and cur_element[0] == "EVAL"
            ):  
                expr = cur_element[1]  
                # replaces each group of variable of the form ${...} to fill_dict["..."]
                expr = re.sub(r"\${([^{}]+)}", r'fill_dict["\1"]', expr)
                val = eval(expr)
                write_dict(conf_exp, cur_key, val)
                fill_dict[cur_key[-1]] = val
            elif isinstance(cur_element, str) and cur_element== "CURRENT_TIME":
                val = time.strftime("%Y%m%d%H%M%S")
                write_dict(conf_exp, cur_key, val)
                fill_dict[cur_key[-1]] = val
            # format: ["BINARY_SEARCH", LB, UB, EPS, FILE, KEYS, EXPR]
            # find the minimal value such that EXPR > 0
            elif (
                isinstance(cur_element, list)
                and isinstance(cur_element[0], str)
                and cur_element[0] == "BINARY_SEARCH"
            ):  
                assert(len(cur_element) == 7)
                lb = cur_element[1] 
                ub = cur_element[2]
                eps = cur_element[3]
                filename = cur_element[4]
                keys = cur_element[5]
                expr = cur_element[6]
                if not isinstance(eps, numbers.Number):
                    eps = 1e-6 # default epsilon
                
                if not isinstance(keys, list):
                    # everything is key
                    keys = [x[-1] for x in key_list]
                while ub - lb > eps:
                    mid = (ub + lb) / 2
                    write_dict(conf_exp, cur_key, mid)
                    fill_dict[cur_key[-1]] = mid
                    check_state = {}
                    for key in keys:
                        check_state[key] = []
                    output_format = replace_str(read_dict(conf, filename))
                    expand_recur(conf, conf_exp, key_list, level + 1, fill_dict)
                    # wait for all processes to finish
                    processes.difference_update([p for p in processes if p.poll() is not None])
                    if len(processes) >= 0:
                        os.wait()
                        processes.difference_update([p for p in processes if p.poll() is not None])
                    time.sleep(0.5) # wait for all dumps to be written
                    for f in glob.glob(output_format):
                        with open(f) as fin:
                            output = json.load(fin)
                            for key in keys:
                                check_state[key].append(output[key])
                    expr = re.sub(r"\${([^{}]+)}", r'check_state["\1"]', expr)
                    val = eval(expr)
                    print(val)
                    if val <= 0:  # success
                        lb = mid
                    else:
                        ub = mid
                return
            else:  # do nothing, just update fill_dict
                fill_dict[cur_key[-1]] = cur_element

            level += 1

        # end of recursion
        for key in key_list:
            cur_element = read_dict(conf, key)
            cur_exp_element = read_dict(conf_exp, key)
            if isinstance(cur_element, str):
                new_string = replace_str(cur_element)  
                # so original texts from conf are needed
                write_dict(conf_exp, key, new_string)
            elif isinstance(
                cur_exp_element, list
            ):  # ITER_*** lists in conf are not lists in conf_exp
                for i, x in enumerate(cur_element):  # pattern should keep the same
                    if isinstance(x, str):
                        cur_exp_element[i] = replace_str(x)  # same as above, x is from conf

        # config/auto/auto seems si**y
        par_path = (
            config_folder if config_folder.name == "auto" else config_folder / "auto"
        )
        par_path.mkdir(parents=True, exist_ok=True)
        auto_json_path = str(par_path / "{}_autogen_{}.json".format(config_name, cnt))
        cnt += 1
        with open(auto_json_path, "w") as fout:
            json.dump(conf_exp, fout, indent=4)

        # since the "script" field is not used in the actual program
        # we can update it after the json file is written
        for idx, command in enumerate(conf_exp["script"]):
            if command == "JSON":
                conf_exp["script"][idx] = auto_json_path

        # guard against out of memory
        if psutil.virtual_memory().percent >= max_memory:
            os.wait()
            processes.difference_update([p for p in processes if p.poll() is not None])

        print("running: ", " ".join(conf_exp["script"]))
        processes.add(subprocess.Popen(conf_exp["script"]))

        if len(processes) >= max_processes:
            os.wait()
            processes.difference_update([p for p in processes if p.poll() is not None])
        return

    key_list = []
    expand_keys(config, key_list, [])
    expand_recur(config, copy.deepcopy(config), key_list, 0, {})
    for p in processes:
        p.communicate()


for conf_path in arguments.conf:
    print("Config reading: ", conf_path)
    run_conf(conf_path)
