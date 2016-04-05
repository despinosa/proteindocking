#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .alpslayer import ALPSLayer
from .chromosome import Chromosome
from .definitions.crossover import single_point
from .definitions.fitness import random_fitness
from .definitions.selection import plus
from multiprocessing import Pool

def stop_condition(alpslayer):
    return alpslayer.generation > 250

class ALPS(object):
    def __init__(self):
        Chromosome.setup(random_fitness, {1:'a', 2:'b', 3:'c'})
        ALPSLayer.setup(200, 0.1, 0.8, 5, plus, single_point)
        self.lay0 = ALPSLayer(10)
        self.lay1 = ALPSLayer(20)
        self.lay2 = ALPSLayer(30)
        self.lay3 = ALPSLayer(40)
        self.lay4 = ALPSLayer(50)
        self.lay1.prev_layer = self.lay0
        self.lay2.prev_layer = self.lay1
        self.lay3.prev_layer = self.lay2
        self.lay4.prev_layer = self.lay3
        self.lay0.next_layer = self.lay1
        self.lay1.next_layer = self.lay2
        self.lay2.next_layer = self.lay3
        self.lay3.next_layer = self.lay4

    def run(self):
        while not stop_condition(self.lay0):
            self.lay0.start()
            self.lay1.start()
            self.lay2.start()
            self.lay3.start()
            self.lay4.start()
            self.lay0.join()
            self.lay1.join()
            self.lay2.join()
            self.lay3.join()
            self.lay4.join()
            try:
                print "best: " + self.lay4.population[0]
            except IndexError:
                print "nada"

if __name__ == '__main__':
    a = ALPS()
    a.run()
