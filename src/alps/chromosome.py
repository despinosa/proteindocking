#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scipy as sp
# import numpy as np

sp_rand = sp.random.random
sp_randint = sp.random.randint

class Chromosome(sp.ndarray):
    def __new__(cls, layer, birth, pieces=None):
        # if input_array: arr = np.asarray(input_array).view(cls)
        if pieces is not None:
            arr = sp.concatenate([piece[:] for piece in pieces]).view(cls)
            if arr.size != layer.main.lower.size:
                raise TypeError("longitud incompatible: %d/%d".
                                format(arr.size, layer.main.lower.size))
        else:
            arr = sp_rand((layer.main.lower.size,)).view(cls)
            arr *= layer.main.upper - layer.main.lower
            arr += layer.main.lower
        arr.layer = layer
        arr.birth = birth
        arr.score = arr.fitness()
        return arr

    def __lt__(self, other):
        return self.score < other.score
    def __le__(self, other):
        return self.score <= other.score
    def __gt__(self, other):
        return self.score > other.score
    def __ge__(self, other):
        return self.score >= other.score
    def __nonzero__(self):
        return True

    def __eq__(self, other):
        self.flags.writeable = False
        self_hash = hash(self.data)
        self.flags.writeable = True
        other.flags.writeable = False
        other_hash = hash(other.data)
        other.flags.writeable = True
        return self_hash == other_hash

    def __str__(self):
        return "score={:3f}, birth={}, {}".format(self.score, self.birth,
                                                  super(Chromosome, self).
                                                  __str__())

    def __array_finalize__(self, obj):
        if obj is None: return
        self.birth = getattr(obj, 'birth', None)
        self.score = getattr(obj, 'score', None)
        self.hash = getattr(obj, 'hash', None)

    fitness = lambda self: self.layer.main.fitness(self, self.layer)

    def mutate(self):
        idx = sp_randint(0, self.size)
        self[idx] = self.layer.main.lower[idx] + sp_rand() * (self.layer.main.upper[idx]-
                                                        self.layer.main.lower[idx])
        self.score = self.fitness()
