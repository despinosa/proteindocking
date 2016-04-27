#!/usr/bin/env python
# -*- coding: utf-8 -*-

def plus(layer, new_pop, old_pop=[], **kwargs):
    """"Elitismo para lamdba más mu

    Devuelve los mejores elementos de la unión de las poblaciones dadas
    en tiempo O((lambda+mu)*log(lambda+mu)).

    pop_size ~ tamaño de la población seleccionada
    new_pop ~ población de tamaño lambda de individuos recién generados
    old_pop ~ población de tamaño mu de individuos viejos
    """
    population = new_pop + old_pop
    population.sort()
    del population[layer.main.pop_size:]
    return population


def comma(layer, new_pop, **kwargs):
    """Elitismo para lamdba coma mu

    Devuelve los mejores elementos de la población nueva en tiempo
    O(1).

    pop_size ~ tamaño de la población seleccionada
    new_pop ~ población de tamaño lambda de individuos recién generados
    old_pop ~ población de tamaño mu de individuos viejos
    """
    del new_pop[layer.main.pop_size:]
    return new_pop


def enhanced(layer, new_pop, old_pop=[]):
    del old_pop[layer.main.preserve:]
    create = layer.main.pop_size - len(old_pop)
    del new_pop[create:]
    population = new_pop + old_pop
    population.sort()
    return population
