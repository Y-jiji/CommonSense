# Vector kernel Doro sketch
# "doro" means "mud" in Japanese

import numpy as np
from oniakHash import ohash
import random

# random.seed(3000750715)


class Doro:
    # initialize an empty code
    # By default, counters are int8 types, so values interted must be integers.
    # if count_min is True, then all signs are 1
    def __init__(self, size, k, seeds=None, count_min=False):
        self.array = np.zeros(size, dtype=np.int8)
        # value is stored in a dictionary
        self.value = {}
        self.value_ori = {}
        if seeds is None:
            seeds = [random.randint(0, 0xFFFFFFFF) for _ in range(k)]
        # k hash functions that maps elements into array
        self.hashfunc = [ohash.WYHash(seeds[i], 2 * size) for i in range(k)]

        # summaries of peeling results
        self.num_peels = 0
        self.num_correct_peels = 0
        self.k = k
        self.count_min = count_min

    def reset(self):
        self.array = np.zeros(len(self.array), dtype=np.int8)
        self.value = {}
        self.num_peels = 0
        self.num_correct_peels = 0
        self.encode(self.setD_ori)

    def encode(self, kvpairs: dict | set):
        assert self.value == {}
        self.value = kvpairs
        self.setD_ori = kvpairs.copy()

        if isinstance(kvpairs, set):
            for element in kvpairs:
                for i in range(self.k):
                    index, sign = self.hash(i, element)
                    self.array[index] += sign
        else:
            for element, val in kvpairs.items():
                for i in range(self.k):
                    index, sign = self.hash(i, element)
                    self.array[index] += sign * val

    def hash(self, i, element):
        index = self.hashfunc[i].hash(element)
        sign = 1 if self.count_min or index >= len(self.array) else -1
        index = index % len(self.array)
        return index, sign

    # if delta is 1, all counters will be REDUCED by 1 (if positive sign)
    # or INCREASED by 1 (if negative sign)
    def peel(self, element, delta):
        self.num_peels += 1
        for i in range(self.k):
            index, sign = self.hash(i, element)
            self.array[index] -= sign * delta

        ori_value = self.value.get(element, 0)
        new_value = ori_value - delta
        self.value[element] = new_value
        if abs(new_value) < abs(ori_value):
            self.num_correct_peels += 1

    # senses the signal of an element by inner product
    # performance bottleneck
    def sense(self, element):
        signal = np.zeros(len(self.array), dtype=np.int8)
        for i in range(self.k):
            index, sign = self.hash(i, element)
            signal[index] = sign
        return np.inner(self.array, signal)

    # print result summary
    def show_result(self):
        print("Number of peels: ", self.num_peels)
        print("Number of correct peels: ", self.num_correct_peels)
        print("Number of wrong peels: ", self.num_peels - self.num_correct_peels)
        print("Size of nonzero indices:", len(self.value))
        print("L1 Norm of indices", sum(abs(x[1]) for x in self.value.items()))
