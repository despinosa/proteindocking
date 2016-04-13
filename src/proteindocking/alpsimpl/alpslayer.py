#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bisect import insort
from chromosome import Chromosome
from heapq import nsmallest
from itertools import repeat
from math import ceil
from random import randint, random, sample
from sys import maxint
from threading import Event,  Lock, Thread


class ALPSLayer(Thread):
    def __init__(self, max_age, prev_layer=None, next_layer=None):
        super(ALPSLayer, self).__init__()
        self.name = "{}{}".format(self.__class__.__name__, self.count)
        self.__class__.count += 1
        self.max_age = max_age
        self.prev_layer = prev_layer
        self.next_layer = next_layer
        self.population = []
        self.copied = Event()
        self.copied.clear()
        self.redisted = Event()
        self.redisted.clear()

    @classmethod
    def setup(cls, pop_size, mutate_rate, mating_rate, tourn_size,
              stop_condition, elitism,  crossover, n_parents=2):
        cls.pop_size = pop_size
        cls.stop_condition = classmethod(stop_condition)
        cls.crossover = staticmethod(crossover)
        cls.elitism = staticmethod(elitism)
        cls.tourn_size = tourn_size
        cls.mutate_cycles = int(ceil(mutate_rate * pop_size))
        cls.reprod_cycles = int(ceil(mating_rate * pop_size / 2))
        cls.n_parents = n_parents
        cls.generation = 0
        cls.count = 0

    @classmethod
    def stop_condition(cls):
        raise NotImplementedError

    @staticmethod
    def crossover(*args):
        raise NotImplementedError

    @staticmethod
    def elitism(*args):
        raise NotImplementedError

    n_parents = 2
    generation = 0
    count = 0


    def rand_pop(self):
        del self.population[:]
        for _ in repeat(None, self.pop_size):
            insort(self.population, Chromosome(self.generation))


    def redistribute(self):
        if self.prev_layer is None:
            self.__class__.generation += 1
            if self.generation % self.max_age == 0:
                while len(self.population) > 0:
                    insort(self.next_layer.population, self.population.pop())
                self.rand_pop()
        else:
            i = 0
            while i < len(self.population):
                if self.generation - self.population[i].birth >= self.max_age:
                    insort(self.next_layer.population, self.population[i])
                    self.population.remove(i)                                   
                i += 1
        del self.population[self.pop_size:] # trim


    def iterate(self):
        def reproduce(pool):
            offspring = []
            if len(pool) >= self.tourn_size:
                for _ in repeat(None, self.reprod_cycles):
                    tournament = sample(pool, self.tourn_size)
                    parents = nsmallest(self.n_parents, tournament)
                    for child in self.crossover(*parents):
                        insort(offspring, child) # n
            return offspring

        def mutate(offspring):
            off_size = len(offspring)
            if off_size >= 1:
                for _ in repeat(None, self.mutate_cycles): # mutación
                    mutated = offspring.pop(randint(0, off_size-1))
                    mutated.mutate()
                    insort(offspring, mutated)
            return offspring

        pool = self.population[:]
        if self.prev_layer is not None: # and self.generation > self.min_age:
            pool += self.prev_layer.population[:] #! bloquear
            self.copied.set()
        offspring = reproduce(pool)
        offspring = mutate(offspring)
        if self.next_layer is not None:
            self.next_layer.copied.wait()
            self.next_layer.copied.clear()
        self.population = self.elitism(self.pop_size, offspring,
                                       self.population)
        if self.prev_layer is not None:
            self.prev_layer.redisted.wait()
            self.prev_layer.redisted.clear()
        if self.next_layer is not None:
            self.redistribute()
            self.redisted.set()


    def run(self):
        if self.prev_layer is None:
            self.rand_pop()
        while not self.stop_condition():
            self.iterate()
            if self.next_layer is None:
                if len(self.population)>0:                     
                    print self.population[0]
        self.copied.set()
        self.redisted.set()

