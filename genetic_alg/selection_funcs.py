from qheap import heapify

def plus(pop_size, new_pop, old_pop=None):
    """"elitism over lamdba plus mu"""
    population = new_pop + old_pop
    heapify(population)
    return population[:pop_size]

def comma(pop_size, new_pop, old_pop=None):
    """elitism over lamdba comma mu"""
    return new_pop[:pop_size]
