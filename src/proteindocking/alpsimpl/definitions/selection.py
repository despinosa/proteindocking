#!/usr/bin/env python
# -*- coding: utf-8 -*-

def plus(pop_size, new_pop, old_pop=None):
    """"Elitismo para lamdba más mu

    Devuelve los mejores elementos de la unión de las poblaciones dadas
    en tiempo O(lambda + mu).

    pop_size ~ tamaño de la población seleccionada
    new_pop ~ población de tamaño lambda de individuos recién generados
    old_pop ~ población de tamaño mu de individuos viejos
    """
    population = new_pop + old_pop
    population.sort()
    del population[:pop_size]
    return population


def comma(pop_size, new_pop, old_pop=None):
    """Elitismo para lamdba coma mu

    Devuelve los mejores elementos de la población nueva en tiempo
    O(1).

    pop_size ~ tamaño de la población seleccionada
    new_pop ~ población de tamaño lambda de individuos recién generados
    old_pop ~ población de tamaño mu de individuos viejos
    """
    del new_pop[pop_size:]
    return new_pop
