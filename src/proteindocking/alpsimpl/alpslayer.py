#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bisect import insort
from chromosome import Chromosome
from heapq import nsmallest
from itertools import repeat
from math import ceil
from random import randint, random, sample
from sys import maxint
from threading import Event, Thread


class ALPSLayer(Thread):
    def __init__(self, max_age, prev_layer=None, next_layer=None):
        super(ALPSLayer, self).__init__()
        self.max_age = max_age
        if prev_layer is not None:
            self.min_age = prev_layer.max_age
        else:
            self.min_age = 0
        self.first = None
        self.age_gap = self.max_age - self.min_age
        self.prev_layer = prev_layer
        self.next_layer = next_layer
        self.population = []
        self.ready = Event()
        self.ready.clear()

    @classmethod
    def setup(cls, pop_size, mutate_rate, mating_rate, tourn_size, elitism,
              crossover, n_parents=2):
        cls.pop_size = pop_size
        cls.crossover = staticmethod(crossover)
        cls.elitism = staticmethod(elitism)
        cls.tourn_size = tourn_size
        cls.mutate_cycles = int(ceil(mutate_rate * pop_size))
        cls.reprod_cycles = int(ceil(mating_rate * pop_size / 2))
        cls.n_parents = n_parents

    @staticmethod
    def crossover(*args):
        raise NotImplementedError

    @staticmethod
    def elitism(*args):
        raise NotImplementedError

    generation = 1
    n_parents = 2

    def rand_pop(self):
        self.population = []
        for _ in repeat(self.pop_size):
            insort(self.population, Chromosome())

    def iterate(self):
        def reproduce(source):
            offspring = []
            source_size = len(source)
            if source_size > self.n_parents:
                for _ in repeat(None, self.reprod_cycles):
                    tournament = sample(source, min(self.tourn_size,
                                                    source_size))
                    parents = nsmallest(self.n_parents, tournament)
                    for child in self.crossover(*parents):
                        insort(offspring, child)
            return offspring

        def mutate(offspring):
            off_size = len(offspring)
            if off_size > 0:
                for _ in repeat(None, self.mutate_cycles): # mutación
                    mutated = offspring.pop(randint(0, off_size-1))
                    mutated.mutate()
                    insort(offspring, mutated)
            return offspring

        def redistribute():
            self.next_layer.population += self.population
            self.next_layer.population.sort()
            if self.prev_layer is None:
                self.rand_pop()
            else:
                self.population = []

        if self.min_age is None and self.prev_layer is not None:
            self.min_age = (self.prev_layer.max_age if self.prev_layer is None
                                                    else 0)
        if self.age_gap is None:
            self.age_gap = (self.max_age - self.min_age if self.next_layer is
                                                        not None else maxint)
        source = self.population
        if self.prev_layer is not None and self.generation > self.min_age:
            if self.first == None:
                self.first = self.generation
            source += self.prev_layer.population #! bloquear
        self.ready.set()
        offspring = reproduce(source)
        offspring = mutate(offspring)
        if self.next_layer is not None:
            self.next_layer.ready.wait()
        self.population = self.elitism(self.pop_size, offspring,
                                       self.population)
        if (self.generation % self.age_gap == self.first and
                self.next_layer is not None):
            redistribute()
        del self.population[self.pop_size:] # trim


    run = iterate


    def join(self):
        super(ALPSLayer, self).join()
        super(ALPSLayer, self).__init__()
        self.ready.clear()
