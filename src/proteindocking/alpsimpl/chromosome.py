#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scipy as sp
# import numpy as np

sp_rand = sp.random.random
sp_randint = sp.random.randint
sp_rand_angles = sp.random.uniform

class Chromosome(sp.ndarray):
    @classmethod
    def setup(cls, fitness, encoding):
        cls.fitness = fitness
        cls.encoding = encoding
        cls.length = len(encoding)

    def __new__(cls, birth, pieces=None):
        # if input_array: arr = np.asarray(input_array).view(cls)
        if pieces is not None:
            arr = sp.concatenate([piece[:] for piece in pieces]).view(cls)
            if arr.length != cls.length:
                raise TypeError("longitud incompatible: %d/%d".
                                format(arr.length, cls.length))
        else:
            arr = sp_rand((cls.length,)).view(cls)
            #arr = sp_randint(10,size=2).view(cls)
            #arr = sp_rand_angles(0,360,cls.length).view(cls)
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

    def fitness(self):
        raise NotImplementedError

    encoding = None
    length = None

    def mutate(self):
        self[sp_randint(0, Chromosome.length)] = sp_rand()
        # self[sp_randint(0, Chromosome.length)] = sp_rand_angles(0,360)
        self.score = self.fitness()
