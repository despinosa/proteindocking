#!/usr/bin/env python
# -*- coding: utf-8 -*-

from alpslayer import ALPSLayer
from chromosome import Chromosome
from definitions.crossover import single_point
from definitions.fitness import rastrigin
from definitions.selection import plus

GEN_LIMIT = 500

def stop_condition(cls):
    return cls.generation >= GEN_LIMIT

class ALPS(object):
    def __init__(self):
        Chromosome.setup(rastrigin, {1:'a', 2:'b', 3:'b', 4:'d', 5:'e'})
        ALPSLayer.setup(50, 0.1, 0.8, 5, stop_condition, plus, single_point)
        self.layers = [ALPSLayer(10*i) for i in xrange(1, 11)]
        for i in xrange(1, len(self.layers)):
            self.layers[i].prev_layer = self.layers[i-1]
        for i in xrange(len(self.layers)-1):
            self.layers[i].next_layer = self.layers[i+1]

    def run(self):
        for lay in self.layers:
            lay.start()
        for lay in self.layers:
            lay.join()

if __name__ == '__main__':
    a = ALPS()
    a.run()
