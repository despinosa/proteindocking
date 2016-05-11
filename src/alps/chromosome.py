#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from threading import Lock
# import numpy as np

sp_rand = np.random.random
sp_randint = np.random.randint

class Chromosome(np.ndarray):
    def __new__(cls, main, birth, pieces=None):
        # if input_array: arr = np.asarray(input_array).view(cls)
        if pieces is not None:
            arr = np.concatenate([piece[:] for piece in pieces]).view(cls)
            if arr.size != main.lower.size:
                raise TypeError("longitud incompatible: %d/%d".
                                format(arr.size, main.lower.size))
        else:
            arr = sp_rand((main.lower.size,)).view(cls)
            arr *= main.span
            arr += main.lower
        arr.main = main
        arr.birth = birth
        arr.score = None
        arr.flags.writeable = False
        arr.hash = hash(arr.data)
        arr.lock = Lock()
        return arr

    def __lt__(self, other):
        if self.score is None: self.eval_fitness()
        if other.score is None: other.eval_fitness()
        return self.score < other.score

    def __le__(self, other):
        if self.score is None: self.eval_fitness()
        if other.score is None: other.eval_fitness()
        return self.score <= other.score

    def __gt__(self, other):
        if self.score is None: self.eval_fitness()
        if other.score is None: other.eval_fitness()
        return self.score > other.score

    def __ge__(self, other):
        if self.score is None: self.eval_fitness()
        if other.score is None: other.eval_fitness()
        return self.score >= other.score

    def __nonzero__(self):
        return True

    def __eq__(self, other):
        return self.hash == other.hash

    def __str__(self):
        return "score={0:3f}, birth={1}, {2}".format(self.score, self.birth,
                                                     super(Chromosome, self).
                                                             __str__())

    def __array_finalize__(self, obj):
        if obj is None: return
        self.birth = getattr(obj, 'birth', None)
        self.score = getattr(obj, 'score', None)
        self.hash = getattr(obj, 'hash', None)
        self.lock = getattr(obj, 'lock', None)

    def eval_fitness(self):
        with self.lock:
            self.score = self.main.fitness(self)
        return self.score

    def mutate(self):
        idx = sp_randint(0, self.size)
        with self.lock:
            self.flags.writeable = True
            self[idx] = self.main.lower[idx] + sp_rand()*(self.main.upper[idx]-
                                                          self.main.lower[idx])
            self.flags.writeable = False
            self.hash = hash(self.data)
