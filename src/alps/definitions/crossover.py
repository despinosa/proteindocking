#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ..chromosome import Chromosome
from random import randint

def single_point(alps, parent1, parent2):
    assert parent1.size == parent2.size
    point = randint(0, parent1.size)
    birth = min(parent1.birth, parent2.birth)
    offspring1 = Chromosome(alps, birth, (parent1[:point], parent2[point:]))
    offspring2 = Chromosome(alps, birth, (parent2[:point], parent1[point:]))
    return offspring1, offspring2
