#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ..chromosome import Chromosome
from random import randint

def single_point(parent1, parent2):
    assert parent1.length == parent2.length
    point = randint(0, parent1.length)
    birth = min(parent1.birth, parent2.birth)
    offspring1 = Chromosome(birth, (parent1[:point], parent2[point:]))
    offspring2 = Chromosome(birth, (parent2[:point], parent1[point:]))
    return offspring1, offspring2
