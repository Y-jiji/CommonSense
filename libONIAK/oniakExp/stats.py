import numpy as np

"""
    Result length is nbins.
    Each value corresponding to x is Pr[X <= x].
"""
def cdf(data, limits, nbins):
    data = np.sort(data)
    result = np.empty(nbins, dtype=np.float32)
    xs = np.linspace(limits[0], limits[1], nbins)
    cnt = 0
    idx = 0

    for x in data:
        while x > xs[idx]:
            result[idx] = cnt
            idx += 1
            if idx == nbins:
                break  # break all loops
        else:
            cnt += 1
            continue
        break

    while idx < nbins:
        result[idx] = cnt
        idx += 1
        
    return result / len(data)

def icdf(data, nbins):
    data = np.sort(data)
    ys = np.linspace(0, 0.99999, nbins)
    ys = (ys * len(data)).astype(np.int32)
        
    return data[ys]