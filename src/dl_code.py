import numpy as np

class DLCode:
    # initialize an empty code
    def __init__(self, size):
        self.array = np.zeros(size, dtype=np.int8)
        self.set = set()
        # hash function that maps elements into array
        self.hashfunc = None

        # summaries of peeling results
        self.num_peels = 0
        self.num_correct_peels = 0


    # D = A \ B
    def encode(self, setD : set):
        assert(self.set == set())
        self.set = setD

        # TODO: encode the set into array
        pass
    
    def peel(self, element, direction):
        pass
    
    # pring result summary
    def show_result(self):
        pass

    
    