#!/usr/bin/env python
# -*- coding: utf-8 -*-

from alpslayer import ALPSLayer
from chromosome import Chromosome
from definitions.crossover import single_point
from definitions.fitness import random_fitness
from definitions.selection import plus

def stop_condition(cls):
    return cls.generation >= 1000

class ALPS(object):
    def __init__(self):
        Chromosome.setup(random_fitness, {1:'a', 2:'b', 3:'c'})
        ALPSLayer.setup(50, 0.1, 0.8, 5, stop_condition, plus, single_point)
        self.layers = [ALPSLayer(2**i) for i in xrange(1, 4)]
        for i in xrange(1, len(self.layers)):
            self.layers[i].prev_layer = self.layers[i-1]
        for i in xrange(len(self.layers)-1):
            self.layers[i].next_layer = self.layers[i+1]

    def run(self):
        for lay in self.layers:
            lay.start()
        # try:
        #     print i, lay.population[0]
        # except IndexError:
        #     pass
        for lay in self.layers:
            lay.join()
            print 'termina {}'.format(lay.name)
            print 'estado {} {}\n'.format(lay.copied.is_set(),
                                          lay.redisted.is_set())

if __name__ == '__main__':
    a = ALPS()
    a.run()
