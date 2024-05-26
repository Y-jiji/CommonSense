from oniakIO import odats, ojson
import sys, random, time, json
from cdma_src.dl_code import *
from cdma_src.decode import *
import os

"""
        Terminologies: B is the set of message sender.
        A is the set of message receiver.
        A recovers A \ B from the reconciliation process.
        C is the set of recovered elements from the code, 
        which may contain correct or wrong recoveries.
        D is the set of false negatives.
        E is the set of false positives.
"""


config = ojson.load_config(sys.argv[1])
universe = config["universe"]
a_size = config["A size"]
b_size = config["B size"]
a_minus_b_size = config.get("A minus B size", a_size - b_size)
k = config["k"]
d = config["d"]
t0 = config.get("t0", k // 2)
tk = config.get("tk", 20)
save_path = config["result filename"]
max_rounds = config.get("max rounds", 100)
countmin = config.get("countmin", False)

if os.path.exists(save_path):
    print("File already exists")
    sys.exit(1)

listU = list(range(universe))
random.shuffle(listU)
setA = set(listU[:a_size])
b_end = a_minus_b_size + b_size
setB = set(listU[a_minus_b_size:b_end])
assert b_end >= a_size
b_minus_a_size = b_end - a_size 

code = DLCode(d, k, count_min=countmin)
code.encode(setA.symmetric_difference(setB))
controller = RecentSuccessController(10, 8, 1000)
decoder = DLCDecoder(controller)  # controller is not used for now
stats = {}

start = time.time()
num_rounds = decoder.decode2(code, setA, t0=t0, tk=tk, max_rounds=max_rounds, stats=stats)
end = time.time()

result = {}
result["time"] = end - start
result["success"] = len(code.setD) + len(code.setE) == 0
result["num peels"] = code.num_peels
result["num correct peels"] = code.num_correct_peels
result["num rounds"] = num_rounds
result["final D size"] = len(code.setD) - b_minus_a_size
result["final E size"] = len(code.setE)

with open(save_path, "w") as f:
    json.dump(result, f, indent=4)
