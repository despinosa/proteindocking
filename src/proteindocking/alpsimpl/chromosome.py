#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scipy as sp

sp_rand = sp.random.random
sp_randint = sp.random.randint

class Chromosome(sp.ndarray):
    @classmethod
    def setup(cls, fitness, encoding):
        cls.fitness = fitness
        cls.encoding = encoding
        cls.length = len(encoding)

    def __new__(cls, age=0, pieces=None):
        # if input_array: arr = np.asarray(input_array).view(cls)
        if pieces is not None:
            arr = sp.concatenate([piece[:] for piece in pieces]).view(cls)
            if arr.length != cls.length:
                raise TypeError("longitud incompatible: %d/%d".
                                format(arr.length, cls.length))
        else:
            arr = sp_rand((cls.length,)).view(cls)
        arr.age = age
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
        return "{:5f}: {}".format(self.score, super(Chromosome, self).__str__())

    def __array_finalize__(self, obj):
        if obj is None: return
        self.age = getattr(obj, 'age', None)
        self.score = getattr(obj, 'score', None)

    def fitness(self):
        raise NotImplementedError

    encoding = None
    length = None
    def increment_age(self):
        self.age += 1
    def mutate(self):
        self[sp_randint(0, Chromosome.length)] = sp_rand()
        self.score = self.fitness()
