#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scipy as sp
# import numpy as np

sp_rand = sp.random.random
sp_randint = sp.random.randint

class Chromosome(sp.ndarray):
    def __new__(cls, main, birth, pieces=None):
        # if input_array: arr = np.asarray(input_array).view(cls)
        if pieces is not None:
            arr = sp.concatenate([piece[:] for piece in pieces]).view(cls)
            if arr.size != main.lower.size:
                raise TypeError("longitud incompatible: %d/%d".
                                format(arr.size, main.lower.size))
        else:
            arr = sp_rand((main.lower.size,)).view(cls)
            arr *= main.upper - main.lower
            arr += main.lower
        arr.main = main
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

    def __str__(self):
        return "score={:3f}, birth={}, {}".format(self.score, self.birth,
                                                  super(Chromosome, self).
                                                  __str__())

    def __array_finalize__(self, obj):
        if obj is None: return
        self.birth = getattr(obj, 'birth', None)
        self.score = getattr(obj, 'score', None)

    fitness = lambda self: self.main.fitness(self)

    def mutate(self):
        idx = sp_randint(0, self.size)
        self[idx] = self.main.lower[idx] + sp_rand() * (self.main.upper[idx]-
                                                        self.main.lower[idx])
        self.score = self.fitness()
