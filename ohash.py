import wyhash, sys

class WYHash:
    def __init__(self, seed, mod):
        self.seed = seed
        self.sec = wyhash.make_secret(seed)
        self.mod = mod

    # assuming key is int32
    def hash(self, key: int):
        binary_key = key.to_bytes(4, sys.byteorder)
        return wyhash.hash(binary_key, self.seed, self.sec) % self.mod