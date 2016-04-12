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
        self.age_gap = self.max_age - self.min_age
        self.prev_layer = prev_layer
        self.next_layer = next_layer
        self.population = []
        self.ready = Event()
        self.ready.clear()

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


    def rand_pop(self):
        self.population = []
        for _ in repeat(None, self.pop_size):
            insort(self.population, Chromosome(self.generation))


    def distribute(self):
        if self.prev_layer is None:
            self.__class__.generation += 1
            if self.generation % self.max_age == 0:
                for i, individual in enumerate(self.population):
                    insort(self.next_layer.population, self.population.pop(i))
                self.rand_pop()
        else:
            for i, individual in enumerate(self.population):
                if self.generation - individual.birth >= self.max_age:
                    insort(self.next_layer.population, self.population.pop(i))


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

        del self.population[self.pop_size:] # trim
        pool = self.population
        if self.prev_layer is not None: # and self.generation > self.min_age:
            pool += self.prev_layer.population[:] #! bloquear
            self.prev_layer.ready.wait()
        self.ready.set()
        offspring = reproduce(pool)
        offspring = mutate(offspring)
        if self.next_layer is not None:
            self.next_layer.ready.wait()
        self.population = self.elitism(self.pop_size, offspring,
                                       self.population)
        if self.next_layer is not None:
            self.distribute()
        self.ready.clear()


    def run(self):
        if self.prev_layer is None:
            self.distribute()
        while not self.stop_condition():
            self.iterate()


    # def join(self):
    #     super(ALPSLayer, self).join()
    #     super(ALPSLayer, self).__init__()
    #     print self.population
