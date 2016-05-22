    #!/usr/bin/env python
# -*- coding: utf-8 -*-

from abc import ABCMeta, abstractmethod
from alpslayer import ALPSLayer
from chromosome import Chromosome
from definitions.agingscheme import linear
from definitions.crossover import single_point
from math import log

class ALPSException(Exception): pass

class ALPS(object):
    __metaclass__ = ABCMeta

    instances = 0

    def __init__(self):
        super(ALPS, self).__init__()
        self.__class__.instances += 1

    def setup(self, pop_size, mutate_rate, mating_rate, stop_condition,
              elitism, aging_scheme, crossover=single_point,
              max_generations=50, n_layers=5, n_parents=2):
        self.name = '{0}{1}'.format(self.__class__.__name__.lower(),
                                    self.instances)
        self.pop_size = pop_size
        self.stop_condition = lambda *a, **kw: stop_condition(self, *a, **kw)
        self.crossover = lambda *a, **kw: crossover(self, *a, **kw)
        self.elitism = lambda *a, **kw: elitism(self, *a, **kw)
        self.mutate_cycles = int(mutate_rate * pop_size)
        self.preserve = int((1-mating_rate) * pop_size)
        self.reprod_cycles = int(mating_rate * pop_size / 2)
        self.n_parents = n_parents
        self.max_generations = max_generations
        self.remaining_gens = self.max_generations
        self.generation = 0
        age_limits = aging_scheme()
        self.layers = [ALPSLayer(self, i, age_limits.next())
                       for i in xrange(n_layers)]
        for i in xrange(1, len(self.layers)):
            self.layers[i].prev_layer = self.layers[i-1]
        for i in xrange(len(self.layers)-1):
            self.layers[i].next_layer = self.layers[i+1]
        self.best = Chromosome(self, 0)
        self.best.eval_fitness()

    @abstractmethod
    def fitness(self, arr):
        pass

    def start_layers(self):
        for layer in self.layers:
            layer.start()

    def join_layers(self):
        for layer in self.layers:
            layer.join()


if __name__ == '__main__':
    a = ALPS()
    a.run()
