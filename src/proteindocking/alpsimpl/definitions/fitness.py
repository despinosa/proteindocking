#!/usr/bin/env python
# -*- coding: utf-8 -*-

from random import random
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
def levi_function(array):
	x,y = decode(array)
	return ((math.sin(3*math.pi*x))**2) + ((x-1)**2) * (1+ (math.sin(3*math.pi*y))**2) + ((y-1)**2) * (1+ (math.sin(2*math.pi*y))**2)

def easom_function(array):
	x = array[0]
	y = array[1]
	return (-1)*math.cos(x)*math.cos(y)*math.exp(-(((x-math.pi)**2) + ((y-math.pi)**2)))
	
def random_fitness(array):
    return random()
