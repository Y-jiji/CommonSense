from random import *
import numpy as np
from oniakIO import odats, ojson
import time, sys, json
from cdma_src.dl_code import *
from cdma_src.compressive_sensing import compressive_solve
from scipy.sparse import csc_array

config = ojson.load_config(sys.argv[1])
universe = config["universe"]
a_size = config["A size"]
b_size = config["B size"]
k = config["k"]
d = config["d"]
# an item is considered to be in A \ B if its value is in [thre_low, thre_up]
thre_low = config.get("lower threshold", 0.9)
thre_up = config.get("upper threshold", 1.1)
set_filename = config.get("set filename", None)
save_path = config["result filename"]

listU = list(range(universe))
random.shuffle(listU)
setA = set(listU[:a_size])
setB = set(listU[:b_size])
sets = np.array(listU[:a_size], dtype=np.int32)
if set_filename is not None:
    odats.write_file(set_filename, sets)

seeds = [random.randint(0, 0xFFFFFFFF) for _ in range(k)]
code = DLCode(d, k)
code.encode(setA - setB)
y = code.array.astype(np.float32)
hashfuncs = code.hashfunc

a_data = np.ones((k * a_size), dtype=np.float32)
a_colidx = np.zeros_like(a_data)
a_rowidx = np.zeros_like(a_data)
cnt = 0
for j, a in enumerate(listU[: a_size]):
    for i in range(k):
        index, sign = code.hash(i, a)
        a_colidx[cnt] = j
        a_rowidx[cnt] = index
        a_data[cnt] = sign
        cnt += 1
A = csc_array((a_data, (a_rowidx, a_colidx)), shape=(d, a_size))

start = time.time()
x = compressive_solve(A, y)
end = time.time()

setC = set([listU[xi] for xi, xd in enumerate(x) if xd > thre_low and xd < thre_up])
setGT = setA - setB
setD = setGT - setC
setE = setC - setGT

result = {}
result["time"] = end - start
result["success"] = (len(setD) + len(setE) == 0)
result["final D size"] = len(setD)
result["final E size"] = len(setE)

x = x[b_size:]
result["numerical error"] = np.linalg.norm(x - np.ones_like(x), 1)

with open(save_path, "w") as f:
    json.dump(result, f, indent=4)