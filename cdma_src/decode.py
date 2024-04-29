from cdma_src.dl_code import DLCode

PeelingDirection = DLCode.PeelingDirection

def reverse_direction(direction : PeelingDirection):
    if direction == PeelingDirection.forward:
        return PeelingDirection.backward
    elif direction == PeelingDirection.backward:
        return PeelingDirection.forward

""" This controller controls the direction of the peeling process."""
class PeelingController:
        def __init__(self):
            self.peeling_direction = PeelingDirection.forward
            pass
        
        def update(self, args):
            # TODO: update the peeling direction based on new knowledge
            pass

        # reset all internal states to decode a new instance of code
        def reset(self):
            pass

# controls the peeling direction by recent numbers or successes
# needs complete knowledge of encoded set (not available in real use)
class RecentSuccessController(PeelingController):
      def __init__(self, window, threshold, max_reverses):
          self.window = window
          self.threshold = threshold
          self.max_reverses = max_reverses

          self.peeling_direction = PeelingDirection.forward
          self.successes = []
          self.last_count = 0
          self.num_reverses = 0
      
      def update(self, code : DLCode):
          success = code.num_correct_peels - self.last_count
          self.last_count = code.num_correct_peels

          self.successes.append(success)
          if len(self.successes) > self.window:
              self.successes.pop(0)
          if len(self.successes) == self.window and sum(self.successes) < self.threshold:
              self.peeling_direction = reverse_direction(self.peeling_direction)
              self.num_reverses += 1
              self.successes = []
              if self.num_reverses >= self.max_reverses:
                  self.peeling_direction = PeelingDirection.end
          if len(code.setD) + len(code.setE) == 0:
              self.peeling_direction = PeelingDirection.end

      def reset(self):
          self.last_count = 0
          self.successes = []
          self.peeling_direction = PeelingDirection.forward

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
            
            self.controller.update(self.code)
            self.code.show_result()
                 
    def forward_peel(self):
        signals = [(self.code.sense(element), element) for element in self.setA]
        signals.sort(reverse=True)
        for _, element in signals:
            self.code.peel(element, PeelingDirection.forward)
            self.controller.update(self.code)
            self.setC.add(element)

            if self.controller.peeling_direction != PeelingDirection.forward:
                break
    
    def backward_peel(self):
        signals = [(self.code.sense(element), element) for element in self.setC]
        signals.sort()
        for _, element in signals:
            self.code.peel(element, PeelingDirection.backward)
            self.controller.update(self.code)
            self.setC.remove(element)

            if self.controller.peeling_direction != PeelingDirection.backward:
                break
    