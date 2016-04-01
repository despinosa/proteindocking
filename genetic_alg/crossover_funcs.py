from random import randint

def single_point(parent1, parent2):
    point = randint(parent1.length)
    offspring1 = parent1[:point] + parent2[point:]
    offspring2 = parent2[:point] + parent1[point:]
    return offspring1, offspring2
