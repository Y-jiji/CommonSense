import numpy as np
import ohash
import random
from enum import Enum

random.seed(3000750715)

# forward means adding elements in A to C
# backward means removing elements from C
# end means stopping the decoding


class DLCode:
    class PeelingDirection(Enum):
        forward = "forward"
        backward = "backward"
        end = "end"

    # initialize an empty code
    def __init__(self, size, k):
        self.array = np.zeros(size, dtype=np.int8)
        # code in array = encode(setD) - encode(setE)
        self.setD = set()
        self.setE = set()
        # k hash functions that maps elements into array
        self.hashfunc = [ohash.WYHash(random.randint(0, 0xffffffff), 2 * size) for _ in range(k)]

        # summaries of peeling results
        self.num_peels = 0
        self.num_correct_peels = 0
        self.k = k

    # D = A \ B
    def encode(self, setD: set):
        assert self.setD == set()
        self.setD = setD

        for element in setD:
            for i in range(self.k):
                index, sign = self.hash(i, element)
                self.array[index] += sign

    def hash(self, i, element):
        index = self.hashfunc[i].hash(element)
        sign = 1 if index >= len(self.array) else -1
        index = index % len(self.array)
        return index, sign

    def peel(self, element, direction: PeelingDirection):
        self.num_peels += 1
        if direction == self.PeelingDirection.forward:
            peeling_sign = 1
            correct_set = self.setD
            wrong_set = self.setE
        else:
            peeling_sign = -1
            correct_set = self.setE
            wrong_set = self.setD

        for i in range(self.k):
            index, sign = self.hash(i, element)
            self.array[index] -= sign * peeling_sign

        if element in correct_set:
            self.num_correct_peels += 1
            correct_set.remove(element)
        else:
            wrong_set.add(element)

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
        print("Number of elements in D: ", len(self.setD))
        print("Number of elements in E: ", len(self.setE))
