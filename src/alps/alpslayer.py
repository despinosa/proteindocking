#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bisect import insort
from chromosome import Chromosome
from heapq import nsmallest
from itertools import repeat
from math import exp, log
from random import randint, random, sample
from threading import Event, Lock, Thread, current_thread

class ALPSLayer(object): # Thread):
    def __init__(self, main, i, max_age, prev_layer=None, next_layer=None):
        super(ALPSLayer, self).__init__()
        self.main = main
        self.name = "{0}layer{1}".format(main.name, i)
        self.max_age = max_age
        self.prev_layer = prev_layer
        self.next_layer = next_layer
        u_ln_u = self.main.pop_size * log(self.main.pop_size)
        takeover_inv = (float(self.max_age) / self.main.max_generations)
        self.tourn_size = int(round(self.main.n_parents * u_ln_u**takeover_inv))
        self.population = []

    def rand_pop(self):
        del self.population[:]
        for _ in repeat(None, self.main.pop_size):
            insort(self.population, Chromosome(self.main, self.main.generation))
    

    def redistribute(self):
        if self.prev_layer is None:
            self.main.log.write('{}\t{}\n'.format(self.main.generation,
                                                  self.main.best.score))
            self.main.generation += 1
            self.main.remaining_gens -= 1
            if self.main.generation % self.max_age == 0:
                while len(self.population) > 0:
                    insort(self.next_layer.population, self.population.pop())
                self.rand_pop()
        else:
            i = 0
            while i < len(self.population):
                if (self.main.generation - self.population[i].birth >=
                        self.max_age):
                    insort(self.next_layer.population, self.population.pop(i))
                else:
                    i += 1
            del self.population[self.main.pop_size:] # trim

    def iterate(self):
        def reproduce(pool):
            offspring = []
            pool_size = len(pool)
            # tourn_size = 5
            # if pool_size < tourn_size: return offspring
            if pool_size < self.main.n_parents: return offspring
            for _ in repeat(None, self.main.reprod_cycles):
                tournament = sample(pool, min(self.tourn_size, pool_size))
                parents = nsmallest(self.main.n_parents, tournament)
                for child in self.main.crossover(*parents):
                    offspring.append(child) # n
            return offspring

        def mutate(offspring):
            off_size = len(offspring)
            if off_size >= 1:
                for _ in repeat(None, self.main.mutate_cycles): # mutación
                    mutated = offspring.pop(randint(0, off_size-1))
                    mutated.mutate()
                    offspring.append(mutated)
            return offspring

        pool = self.population[:]
        if self.prev_layer is not None:
            if self.main.generation > self.prev_layer.max_age:
                pool += self.prev_layer.population[:] #! bloquear
        offspring = reproduce(pool)
        offspring = mutate(offspring)
        offspring.sort()
        self.population = self.main.elitism(offspring, self.population)
        if self.next_layer is not None:
            self.redistribute()
        if len(self.population) > 0:
            self.main.best = min(self.main.best, self.population[0])
