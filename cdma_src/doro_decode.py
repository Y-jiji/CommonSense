from cdma_src.doro_sketch import Doro


class DoroDecoder:
    def __init__(self):
        self.code = None
        # this is the "suspect" index set in which elements can be nonzero
        self.setA = set()
        # the running result of the decoder
        self.result = {}

    """ current element value + delta = new element value"""

    def get_delta(self, delta, value, value_range):
        new_value = value + delta
        if isinstance(value_range, tuple):
            new_value = max(
                value_range[0], new_value
            )  # if new value is smaller than lower bound, set to lower bound
            new_value = min(
                value_range[1], new_value
            )  # if new value is larger than upper bound, set to upper bound
            new_value = round(new_value)
        # if value_range is a set, set value to the closest one in the set
        elif isinstance(value_range, set):
            distances = [(abs(new_value - d), d) for d in value_range]
            new_value = min(distances)[1]
            if new_value < -2:
                print(distances)
        else:
            raise ValueError
        
        result = new_value - value
        # if result <= -2:
        #         print(distances, new_value, value)
        return new_value - value

    def decode(
        self,
        code: Doro,
        setA: set,
        tk,
        t0=None,
        value_range=None,
        ta=5,
        max_rounds=100,
        verbose=False,
        stats=None,
    ):
        if t0 is None:
            t0 = code.k // 2 + 1
        self.code = code
        self.setA = setA
        thrashing = {}

        for rnd in range(max_rounds):
            # sensing stage
            signals = []
            for element in self.setA:
                power = self.code.sense(element)
                cur_element_value = self.result.get(element, 0)
                delta = self.get_delta(power / code.k, cur_element_value, value_range)
                if abs(delta) > 1e-6:
                    signals.append((thrashing.get(element, 0), power, delta, cur_element_value, element))

            # sort from strong signals to weak ones by absolute value
            signals.sort(key=lambda x: (-abs(x[1]), x[0], x[-1]))
            signals = [
                (thrash, power, delta, cur_value, element)
                for thrash, power, delta, cur_value, element in signals
                if abs(power) >= t0
            ]
            siglen = len(signals)
            if siglen > tk:
                signals = signals[:tk]

            if len(signals) == 0:
                break
            min_thrash = min([thrash for thrash, _, _, _, _ in signals])
            signals2 = [x for x in signals if abs(x[1]) > t0 and x[0] < min_thrash + 2]
            if len(signals2) > 0:
                signals = signals2

            i = 0
            finished = False
            for thrash, power, delta, cur_element_value, element in signals:
                if verbose:
                    print(thrash, power, element, cur_element_value, delta)
                self.code.peel(element, delta)
                self.result[element] = cur_element_value + delta
                i += 1

                thrashing[element] = thrash + 1
                if thrashing[element] > ta : # and siglen > len(code.array) // 16:
                    finished = True
                    break
            i = min(i, len(signals))
            signals = signals[i:]
            if finished:
                break

            if verbose:
                print("Round: ", rnd)
                self.code.show_result()

        if stats is not None:
            stats["signals"] = signals
        return rnd
