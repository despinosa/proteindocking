#!/usr/bin/env python
# -*- coding: utf-8 -*-
import Queue
from bisect import insort
from chromosome import Chromosome
from heapq import nsmallest
from itertools import repeat
from math import ceil, exp
from random import randint, random, sample
from sys import exc_info
from threading import Event, Lock, Thread, current_thread
from traceback import format_exception

class ALPSLayer(Thread):
    def __init__(self, main, i, max_age, prev_layer=None, next_layer=None):
        super(ALPSLayer, self).__init__()
        self.main = main
        self.name = "{0}layer{1}".format(main.name, i)
        self.max_age = max_age
        self.prev_layer = prev_layer
        self.next_layer = next_layer
        self.population = []
        self.copied = Event()
        self.copied.clear()
        self.replaced = Event()
        self.replaced.clear()
        self.ex_queue = Queue.Queue()

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
            pool_size = len(pool)
            if pool_size < self.main.n_parents*2: return offspring
            to_takeover = self.main.generation % self.max_age + 1
            tourn_size = min(max(int(ceil(exp(self.main.ln_sum / to_takeover))),
                                 self.main.n_parents * 2),
                             pool_size)
            for _ in repeat(None, self.main.reprod_cycles):
                tournament = sample(pool, tourn_size)
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
        self.copied.set()
        offspring = reproduce(pool)
        offspring = mutate(offspring)
        offspring.sort()
        if self.next_layer is not None:
            self.next_layer.copied.wait()
            self.next_layer.copied.clear()
        if self.prev_layer is not None:
            self.prev_layer.replaced.wait()
            self.prev_layer.replaced.clear()
        self.population = self.main.elitism(offspring, self.population)
        if self.next_layer is not None:
            self.redistribute()
        try: self.main.best = min(self.main.best, self.population[0])
        except IndexError: pass
        self.replaced.set()

    def run(self):
        try:
            if self.prev_layer is None:
                self.rand_pop()
            while not self.main.stop_condition():
                self.iterate()
            self.copied.set()
            self.replaced.set()
        except Exception as e:
            exc_type, exc_value, exc_traceback = exc_info()
            self.ex_queue.put(repr(format_exception(exc_type, exc_value,exc_traceback)) + repr(e))            
