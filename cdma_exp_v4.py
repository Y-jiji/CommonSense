import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

from oniakIO import odats, ojson
import sys, random, time, json
from cdma_src.doro_sketch import *
from cdma_src.doro_decode import *


"""
        V4 setting: A \ B has value 1, B \ A has value -1, 
        only A \ B is known to the decoder.
"""


config = ojson.load_config(sys.argv[1])
universe = config["universe"]
a_size = config["A size"]
b_size = config["B size"]
a_minus_b_size = config.get("A minus B size", a_size - b_size)
k = config["k"]
d = config["d"]
tk = config.get("tk", 20)
t0 = config.get("t0", None)
save_path = config["result filename"]
max_rounds = config.get("max rounds", 100)
counting = config.get("counting", False)
sg_step = config.get("signal step", 1)
value_range = {0, 1}

if os.path.exists(save_path):
    print("File already exists")
    sys.exit(1)

listU = list(range(universe))
random.shuffle(listU)
setA = set(listU[:a_size])
b_end = a_minus_b_size + b_size
setB = set(listU[a_minus_b_size:b_end])
setA_minus_B = set(listU[:a_minus_b_size])
setB_minus_A = set(listU[a_size:b_end])
assert b_end >= a_size
b_minus_a_size = b_end - a_size

code = Doro(d, k, counting=counting)
gt_vector = {}
for i in range(a_minus_b_size):  # A \ B
    gt_vector[listU[i]] = 1
for i in range(a_size, b_end):  # B \ A
    gt_vector[listU[i]] = -1
code.encode(gt_vector)
decoder = DoroDecoder()
stats = {}

start = time.time()
num_rounds = decoder.decode(
    code,
    setA,  # only set A is known
    tk=tk,
    t0=t0,
    max_rounds=max_rounds,
    stats=stats,
    value_range=value_range,
)
end = time.time()

nonzero_items = set([elem for elem, val in code.value.items() if abs(val) > 1e-6])
nonzero_items = nonzero_items - setB_minus_A

cur_sg = k


result = {}
result["time"] = end - start
result["success"] = len(nonzero_items) == 0
result["num peels"] = code.num_peels
result["num correct peels"] = code.num_correct_peels
result["num rounds"] = num_rounds
result["set D size"] = len(nonzero_items.intersection(setA_minus_B))
result["set E size"] = len(nonzero_items - setA_minus_B)

with open(save_path, "w") as f:
    json.dump(result, f, indent=4)
