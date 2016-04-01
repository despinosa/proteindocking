from . import chromosome.Chromosome
from heapq import heapify, heappush, nsmallest
from itertools import repeat
from model import *
from random import randint, random, sample
from threading import Thread

class ALPSLayer(Thread):
    @classmethod
    def setup(cls, pop_size, crossover, selection, mutate_rate, mating_rate,
              tourn_size):
        cls.pop_size = pop_size
        cls.crossover = crossover
        cls.selection = selection
        cls.tourn_size = tourn_size
        cls.mutate_cycles = mutate_rate * pop_size
        cls.reprod_cycles = mating_rate * pop_size // 2

    def __init__(self, max_age, stop_condition, prev_layer=None,
                 next_layer=None):
        self.max_age = max_age
        self.stop_condition = stop_condition
        self.prev_layer = prev_layer
        self.next_layer = next_layer
        self.generation = 1
        self.randomize()

    def randomize(self):
        assert self.prev_layer is None
        self.population = []
        pushrandom = lambda: heappush(self.population, Chromosome())
        map(pushrandom, repeat(pushrandom, self.pop_size))

    def iterate(self):
        offspring = []
        def reproduce():
            if self.prev_layer is not None:
                tournament = sample(self.population +
                                    self.prev_layer.population,
                                    self.tourn_size)
            else:
                tournament = sample(self.population, self.tourn_size)
            parents = nsmallest(2, tournament)
            pushmany = lambda individual: heappush(offspring, individual)
            map(pushmany, self.crossover(*parents))
        def mutate():
            mutated = offspring.pop(randint(offspring.length))
            mutated.mutate()
            heappush(offspring, mutated)
        def redistribute(i, individual):
            if self.max_age < self.generation - individual.birth:
                self.next_layer.population.heappush(self.population.pop(i))
        map(reproduce, repeat(None, self.reprod_cycles))
        map(mutate, repeat(None, self.mutate_cycles))
        self.population = selection(self.pop_size, offspring, self.population)
        self.generation += 1
        if self.next_layer is not None: # ésta no es ultima capa
            map(redistribute, enumerate(self.population))
        if self.generation % self.max_age == 0:
            if self.prev_layer is None: # ésta es la primer capa
                self.randomize()
        self.population = self.population[:self.pop_size]

