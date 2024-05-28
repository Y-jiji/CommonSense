from cdma_src.doro_sketch import Doro


class DoroDecoder:
    def __init__(self):
        self.code = None
        # this is the "suspect" index set in which elements can be nonzero
        self.setA = set()
        # the running result of the decoder
        self.result = {}

    def decode(
        self,
        code: Doro,
        setA: set,
        t0,
        tk,
        delta_range=None,
        ta=5,
        max_rounds=100,
        verbose=False,
        stats=None,
    ):
        self.code = code
        self.setA = setA
        # terminates decoding if an element is added and removed for multiple times
        thrashing = {}

        for rnd in range(max_rounds):
            signals = [(self.code.sense(element), element) for element in self.setA]

            # sort from strong signals to weak ones by absolute value
            signals.sort(reverse=True, key=lambda x: abs(x[0]))
            # threshold is the tk-th strongest signal
            threshold = signals[tk][0]
            # normally t0 is at least k/2, since otherwise
            # peeling only makes the signal stronger
            if threshold < t0:
                threshold = t0  # lower bound at t0

            finished = True
            for value, element in signals:
                if value < threshold:
                    break
                finished = False
                element_id = abs(element)
                delta = value / code.k
                if delta_range is not None:
                    delta = max(
                        delta_range[0], delta
                    )  # if delta is smaller than lower bound, set to lower bound
                    delta = min(
                        delta_range[1], delta
                    )  # if delta is larger than upper bound, set to upper bound

                self.code.peel(element, delta)
                self.result[element_id] = self.result.get(element_id, 0) + delta

                thrashing[element] = thrashing.get(element, 0) + 1
                if thrashing[element] > ta:
                    finished = True
                    break

            if finished:
                break

            if verbose:
                print("Threshold: ", threshold)
                print("Round: ", rnd)
                self.code.show_result()

        if stats is not None:
            stats["signals"] = signals
        return rnd
