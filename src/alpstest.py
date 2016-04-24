from alps.alps import ALPS
from alps.alpslayer import ALPSLayer
from alps.chromosome import Chromosome
from alps.definitions.crossover import single_point
from alps.definitions.selection import enhanced
from alps.definitions.stopcondition import gen_limit
from random import random
import scipy as sp
import math

limiteA=-20
limiteB=20

def decode(array):
    firstVariableDecoded = limiteA+(limiteB-limiteA)*array[0]
    secondVariableDecoded = limiteA+(limiteB-limiteA)*array[1]
    return firstVariableDecoded,secondVariableDecoded

def f(array):
    x,y = decode(array)
    return -20*math.exp(-0.2*math.sqrt(0.5*(x*x+y*y)))-math.exp(0.5*math.cos(2*math.pi*x)+math.cos(2*math.pi*y)) + math.e + 20

def levi(array):
    x,y = decode(array)
    return ((math.sin(3*math.pi*x))**2) + ((x-1)**2) * (1+ (math.sin(3*math.pi*y))**2) + ((y-1)**2) * (1+ (math.sin(2*math.pi*y))**2)

def easom(array):
    x = array[0]
    y = array[1]
    return (-1)*math.cos(x)*math.cos(y)*math.exp(-(((x-math.pi)**2) + ((y-math.pi)**2)))

def rastrigin(array):
    A = 10
    val = A * array.size
    for x in array:
        val += x**2 - A*math.cos(2*math.pi*x)
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
        self.layers = [ALPSLayer(self, i, 10*i) for i in xrange(1, 11)]
        for i in xrange(1, len(self.layers)):
            self.layers[i].prev_layer = self.layers[i-1]
        for i in xrange(len(self.layers)-1):
            self.layers[i].next_layer = self.layers[i+1]

    def run(self):
        for lay in self.layers:
            lay.start()
        for lay in self.layers:
            lay.join()

    def fitness(self, arr):
        return rastrigin(arr)

if __name__ == '__main__':
    from datetime import datetime
    from sys import argv
    a = ALPSTest(int(argv[1]))
    start = datetime.now()
    a.run()
    print 'time: {}\n'.format(datetime.now()-start)
    print 'best: {}\n\n'.format(a.layers[-1].population[0])
