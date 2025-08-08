from oniakIO import odats, ojson
import sys, random, time, json
from cdma_src.dl_code import *
from cdma_src.decode import *
import os

config = ojson.load_config(sys.argv[1])
universe = config["universe"]
a_size = config["A size"]
b_size = config["B size"]
k = config["k"]
d = config["d"]
t0 = config.get("t0", k//2)
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
setB = set(listU[:b_size])

code = DLCode(d, k, count_min=countmin)
code.encode(setA - setB)
controller = RecentSuccessController(10, 8, 1000)
decoder = DLCDecoder(controller)  # controller is not used for now

start = time.time()
num_rounds = decoder.decode2(code, setA, t0=t0, tk=tk, max_rounds=max_rounds)
end = time.time()

result = {}
result["time"] = end - start
result["success"] = (len(code.setD) + len(code.setE) == 0)
result["num peels"] = code.num_peels
result["num correct peels"] = code.num_correct_peels
result["num rounds"] = num_rounds
result["final D size"] = len(code.setD)
result["final E size"] = len(code.setE)

with open(save_path, "w") as f:
    json.dump(result, f, indent=4)