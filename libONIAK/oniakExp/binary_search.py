import numpy as np

"""
    Compute binary search for func(x) == output on multiple output values.
    Function input must be integers.
    Requires both desired output sequence and func() to be monotinically increasing.
"""
class MultiBinarySearch:
    def __init__(self):
        # probability of having at least T successes in L trials
        self.transform = None
        self.lb_recalls = np.empty(0)
        self.ub_recalls = np.empty(0)
        self.lb_thresholds = np.empty(0)
        self.ub_thresholds = np.empty(0)
        self.count = 0
    
    def calc(self, func, outputs, lb=0, ub=2**31-1):
        out_size = len(outputs)
        self.lb_recalls = np.empty(out_size, func(lb))
        self.ub_recalls = np.empty(out_size, func(ub))
        self.lb_thresholds = np.full(out_size, lb)
        self.ub_thresholds = np.full(out_size, ub)
        self.cnt = 0
        for i, ipt in enumerate(outputs):
            ubi = self.ub_thresholds[i]
            lbi = self.lb_thresholds[i]
            
            while ubi - lbi > 1:
                self.cnt += 1
                mid = (ubi + lbi) // 2
                mid_recall = func(mid)

                j = i
                while j < out_size and mid_recall > outputs[j]:
                    j += 1
                # mid is ub for i<j and lb for i>=j
                for k in range(i, j):
                    if mid_recall <= self.ub_recalls[k]:
                        self.ub_thresholds[k] = mid
                        self.ub_recalls[k] = mid_recall
                for k in range(j, out_size):
                     if mid_recall >= self.lb_recalls[k]:
                        self.lb_thresholds[k] = mid
                        self.lb_recalls[k] = mid_recall
                
                ubi = self.ub_thresholds[i]
                lbi = self.lb_thresholds[i]
        return self.ub_thresholds
