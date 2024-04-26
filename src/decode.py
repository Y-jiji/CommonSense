from dl_code import DLCode
from enum import Enum

# forward means adding elements in A to C
# backward means removing elements from C
# end means stopping the decoding
PeelingDirection = Enum('PeelingDirection', ['forward', 'backward', 'end'])

""" This controller controls the direction of the peeling process."""
class PeelingController:
        def __init__(self, args):
            self.peeling_direction = PeelingDirection.forward
            pass
        
        def update(self, args):
            # TODO: update the peeling direction based on new knowledge
            pass

        # reset all internal states to decode a new instance of code
        def reset(self):
            pass
        

class DLCDecoder:
    def __init__(self, controller : PeelingController):
        self.controller = controller
        self.code = None
        self.setA = set()
        # the set of elements already peeled
        # C is always a subset of A, but not necessarily a subset of B
        self.setC = set()

        # additional data members such as sorted list of signals

    def decode(self, code : DLCode, setA : set):
        self.code = code
        self.setA = setA

        while True:
            if self.controller.peeling_direction == PeelingDirection.forward:
                 self.forward_peel()
            elif self.controller.peeling_direction == PeelingDirection.backward:
                 self.backward_peel()
            elif self.controller.peeling_direction == PeelingDirection.end:
                 break
            
            self.controller.update()
            self.code.show_result()
                 
    def forward_peel(self):
        pass
    
    def backward_peel(self):
        pass