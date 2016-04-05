#!/usr/bin/env python
# -*- coding: utf-8 -*-

from chromosome import Chromosome
from heapq import heapify, heappush, merge, nsmallest
from itertools import repeat
from math import ceil
from random import randint, random, sample
from threading import Event, Thread
        

class ALPSLayer(Thread):
    def __init__(self, max_age, prev_layer=None, next_layer=None):
        super(ALPSLayer, self).__init__()
        self.max_age = max_age
        self.prev_layer = prev_layer
        self.next_layer = next_layer
        self.generation = 1
        self.population = []
        self.ready = Event()

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

    n_parents = 2

    def reset_pop(self):
        self.population = []
        for _ in repeat(self.pop_size):
            heappush(self.population, Chromosome())

    def iterate(self):
        def reproduce(source):
            offspring = []
            source_size = len(source)
            if source_size > 2:
                for _ in repeat(None, self.reprod_cycles):
                    tournament = sample(source, min(self.tourn_size,
                                                    source_size))
                    parents = nsmallest(self.n_parents, tournament)
                    for child in self.crossover(*parents):
                        heappush(offspring, child)
            return offspring

        def mutate(offspring):
            off_size = len(offspring)
            if off_size > 0:
                for _ in repeat(None, self.mutate_cycles): # mutación
                    mutated = offspring.pop(randint(0, off_size-1))
                    mutated.mutate()
                    heappush(offspring, mutated)
            return offspring

        def redistribute():
            if self.prev_layer is None: # primer capa
                if self.generation % self.max_age == 0:
                    self.next_layer.population = list(
                        merge(self.next_layer.population, self.population)
                    ) #! bloquear
                    self.reset_pop()
            elif self.next_layer is not None: # capa intermedia
                elders = []
                for i, individual in enumerate(self.population):
                    if self.max_age < self.generation - individual.birth:
                        heappush(elders, self.population.pop(i))
                self.next_layer.population = list(
                    merge(self.next_layer.population, elders)
                ) #! bloquear

        self.ready.clear()
        source = self.population
        if self.prev_layer is not None:
            source += self.prev_layer.population #! bloquear
        self.ready.set()
        offspring = reproduce(source)
        offspring = mutate(offspring)
        if self.next_layer is not None:
            self.next_layer.ready.wait()
        self.population = self.elitism(self.pop_size, offspring,
                                       self.population)
        redistribute()
        self.population = self.population[:self.pop_size] # trim
        self.generation += 1


    run = iterate


    def join(self):
        super(ALPSLayer, self).join()
        super(ALPSLayer, self).__init__()
