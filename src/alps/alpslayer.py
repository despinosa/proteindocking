#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bisect import insort
from chromosome import Chromosome
from heapq import nsmallest
from itertools import repeat
from random import randint, random, sample
from sys import maxint
from threading import Event, Lock, Thread


class ALPSLayer(Thread):
    def __init__(self, main, i, max_age, prev_layer=None, next_layer=None):
        super(ALPSLayer, self).__init__()
        self.main = main
        self.name = "{}{}".format(self.__class__.__name__, i)
        self.max_age = max_age
        self.prev_layer = prev_layer
        self.next_layer = next_layer
        self.population = []
        self.copied = Event()
        self.copied.clear()
        self.replaced = Event()
        self.replaced.clear()

    def rand_pop(self):
        del self.population[:]
        for _ in repeat(None, self.main.pop_size):
            insort(self.population, Chromosome(self.main, self.main.generation))
    

    def redistribute(self):
        if self.prev_layer is None:
            self.main.generation += 1
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
            if len(pool) >= self.main.tourn_size:
                for _ in repeat(None, self.main.reprod_cycles):
                    tournament = sample(pool, self.main.tourn_size)
                    parents = nsmallest(self.main.n_parents, tournament)
                    for child in self.main.crossover(*parents):
                        insort(offspring, child) # n
            return offspring

        def mutate(offspring):
            off_size = len(offspring)
            if off_size >= 1:
                for _ in repeat(None, self.main.mutate_cycles): # mutación
                    mutated = offspring.pop(randint(0, off_size-1))
                    mutated.mutate()
                    insort(offspring, mutated)
            return offspring

        pool = self.population[:]
        if self.prev_layer is not None:
            if self.main.generation > self.prev_layer.max_age:
                pool += self.prev_layer.population[:] #! bloquear
        self.copied.set()
        offspring = reproduce(pool)
        offspring = mutate(offspring)
        if self.next_layer is not None:
            self.next_layer.copied.wait()
            self.next_layer.copied.clear()
        if self.prev_layer is not None:
            self.prev_layer.replaced.wait()
            self.prev_layer.replaced.clear()
        self.population = self.main.elitism(offspring, self.population)
        if self.next_layer is not None:
            self.redistribute()
        self.replaced.set()


    def run(self):
        if self.prev_layer is None:
            self.rand_pop()
        while not self.main.stop_condition():
            self.iterate()
        self.copied.set()
        self.replaced.set()
