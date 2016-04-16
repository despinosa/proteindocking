#!/usr/bin/env python
# -*- coding: utf-8 -*-

def plus(pop_size, new_pop, old_pop=[], *args):
    """"Elitismo para lamdba más mu

    Devuelve los mejores elementos de la unión de las poblaciones dadas
    en tiempo O((lambda+mu)*log(lambda+mu)).

    pop_size ~ tamaño de la población seleccionada
    new_pop ~ población de tamaño lambda de individuos recién generados
    old_pop ~ población de tamaño mu de individuos viejos
    """
    population = new_pop + old_pop
    population.sort()
    del population[pop_size:]
    return population


def comma(pop_size, new_pop, *args):
    """Elitismo para lamdba coma mu

    Devuelve los mejores elementos de la población nueva en tiempo
    O(1).

    pop_size ~ tamaño de la población seleccionada
    new_pop ~ población de tamaño lambda de individuos recién generados
    old_pop ~ población de tamaño mu de individuos viejos
    """
    del new_pop[pop_size:]
    return new_pop


def enhanced(pop_size, new_pop, old_pop=[], elit_rate=0.2):
    preserve = int(pop_size * elit_rate)
    del old_pop[preserve:]
    create = pop_size - len(old_pop)
    del new_pop[create:]
    population = new_pop + old_pop
    population.sort()
    return population
