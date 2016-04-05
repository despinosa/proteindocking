#!/usr/bin/env python
# -*- coding: utf-8 -*-

from alpslayer import ALPSLayer
from chromosome import Chromosome
from definitions.crossover import single_point
from definitions.fitness import random_fitness
from definitions.selection import plus

def stop_condition(alpslayer):
    return alpslayer.generation > 250

class ALPS(object):
    def __init__(self):
        Chromosome.setup(random_fitness, {1:'a', 2:'b', 3:'c'})
        ALPSLayer.setup(200, 0.1, 0.8, 5, plus, single_point)
        self.layers = [ALPSLayer(i*10) for i in xrange(1, 6)]
        for i, lay in enumerate(self.layers[1:]):
            lay.prev_layer = self.layers[i-1]
        for i, lay in enumerate(self.layers[:-1]):
            lay.next_layer = self.layers[i+1]

    def run(self):
        while not stop_condition(self.layers[0]):
            for lay in self.layers:
                lay.start()
            for lay in self.layers:
                lay.join()
            try:
                print "best: " + self.layers[-1].population[0]
            except IndexError:
                pass # print "nada"

if __name__ == '__main__':
    a = ALPS()
    a.run()
