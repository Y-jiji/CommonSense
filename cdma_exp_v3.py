from oniakIO import odats, ojson
import sys, random, time, json
from cdma_src.doro_sketch import *
from cdma_src.doro_decode import *
import os

"""
        V3 setting: A \ B has value 1, B \ A has value -1, 
        both A \ B and B \ A are known to the decoder.
"""


config = ojson.load_config(sys.argv[1])
universe = config["universe"]
a_size = config["A size"]
b_size = config["B size"]
a_minus_b_size = config.get("A minus B size", a_size - b_size)
k = config["k"]
d = config["d"]
tk = config.get("tk", 20)
save_path = config["result filename"]
max_rounds = config.get("max rounds", 100)
counting = config.get("counting", False)
sg_step = config.get("signal step", 1)
value_range = set([-1, 0, 1])

if os.path.exists(save_path):
    print("File already exists")
    sys.exit(1)

listU = list(range(universe))
random.shuffle(listU)
setA = set(listU[:a_size])
b_end = a_minus_b_size + b_size
setB = set(listU[a_minus_b_size:b_end])
setA_minus_B = set(listU[:a_minus_b_size])
assert b_end >= a_size
b_minus_a_size = b_end - a_size

code = Doro(d, k, counting=counting)
gt_vector = {}
for i in range(a_minus_b_size):
    gt_vector[listU[i]] = 1
for i in range(a_size, b_end):
    gt_vector[listU[i]] = -1
code.encode(gt_vector)
decoder = DoroDecoder()
stats = {}

start = time.time()
num_rounds = decoder.decode(
    code,
    setA.union(setB),
    tk=tk,
    ta=10,
    max_rounds=max_rounds,
    stats=stats,
    value_range=value_range,
    verbose=True,
)
end = time.time()

nonzeroes = [code.nonzero_num(code.value)]
maes = [code.mae(code.value)]

cur_sg = k
if "signals" in stats.keys():
    value = code.value.copy()
    for x in stats["signals"]:
        cur_value = value.get(x[1], 0)
        value[x[1]] = cur_value - decoder.get_delta(
            x[0] / code.k, cur_value, value_range
        )
        if abs(x[0]) < cur_sg:
            cur_sg -= sg_step
            nonzeroes.append(code.nonzero_num(value))
            maes.append(code.mae(value))


result = {}
result["time"] = end - start
result["success"] = nonzeroes[0] == 0
result["num peels"] = code.num_peels
result["num correct peels"] = code.num_correct_peels
result["num rounds"] = num_rounds
result["nonzeroes"] = nonzeroes
result["maes"] = maes

with open(save_path, "w") as f:
    json.dump(result, f, indent=4)
