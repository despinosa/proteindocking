from alps.alps import ALPS
from alps.alpslayer import ALPSLayer
from alps.chromosome import Chromosome
from alps.definitions.crossover import single_point
from alps.definitions.selection import enhanced
from alps.definitions.stopcondition import gen_limit
from random import random
import scipy as sp
from math import e, exp, sin, cos, pi

limiteA=-20
limiteB=20

def decode(array):
    firstVariableDecoded = limiteA+(limiteB-limiteA)*array[0]
    secondVariableDecoded = limiteA+(limiteB-limiteA)*array[1]
    return firstVariableDecoded,secondVariableDecoded

def f(array):
    x,y = decode(array)
    return -20*exp(-0.2*sqrt(0.5*(x*x+y*y)))-exp(0.5*cos(2*pi*x)+cos(2*pi*y)) + e + 20

def levi(array):
    x,y = decode(array)
    return ((sin(3*pi*x))**2) + ((x-1)**2) * (1+ (sin(3*pi*y))**2) + ((y-1)**2) * (1+ (sin(2*pi*y))**2)

def easom(array):
    x = array[0]
    y = array[1]
    return (-1)*cos(x)*cos(y)*exp(-(((x-pi)**2) + ((y-pi)**2)))

def rastrigin(array):
    A = 10
    val = A * array.size
    for x in array:
        val += x**2 - A*cos(2*pi*x)
    return val

def random_fitness(array):
    return random()


class ALPSTest(ALPS):
    """docstring for ALPSTest"""
    def __init__(self, n):
        super(ALPSTest, self).__init__()
        self.lower = sp.full(n, -5.12)
        self.upper = sp.full(n, 5.12)
        self.setup(50, 0.1, 0.8, 5, gen_limit, enhanced, single_point)

    def run(self):
        for lay in self.layers:
            lay.start()
        for lay in self.layers:
            lay.join()

    def fitness(self, arr, layer):
        return rastrigin(arr)

if __name__ == '__main__':
    from datetime import datetime
    from sys import argv
    a = ALPSTest(int(argv[1]))
    start = datetime.now()
    a.run()
    print 'time: {}\n'.format(datetime.now()-start)
    print 'best: {}\n\n'.format(a.layers[-1].population[0])
