import random
import numpy as np

class Functions:

	@staticmethod
	def random_chromosome_numpy(L_INDIVIDUAL):				
		chromosome = [np.float] * L_INDIVIDUAL
		for i in xrange(0,L_INDIVIDUAL):
			chromosome[i] = np.around(np.random.uniform(0,360),decimals=3)
		return chromosome

	@staticmethod
	def random_float_value():
		return np.random.uniform(0,360)

	@staticmethod	
	def create_wheel(layer,previouslayer):
		if(previouslayer != None):
			population_aux = layer.population + previouslayer.population
			n_individuals_aux = layer.n_individuals + previouslayer.n_individuals
		else:
			n_individuals_aux = layer.n_individuals
			population_aux = layer.population
		population_aux.sort()
		Lwheel=n_individuals_aux*10
		maxValue=max(population_aux)
		acc=0
		for p in xrange(n_individuals_aux):
			acc+=maxValue.score - population_aux[p].score
		fraction=[]
		for p in xrange(n_individuals_aux):
			fraction.append( float(maxValue.score - population_aux[p].score)/acc)
			if fraction[-1]<=1.0/Lwheel:
				fraction[-1]=1.0/Lwheel	
		fraction[0]-=(sum(fraction)-1.0)/2
		fraction[1]-=(sum(fraction)-1.0)/2
		wheel=[]
		pc=0

		for f in fraction:
			Np=int(f*Lwheel)
			for i in xrange(Np):
				wheel.append(pc)
			pc+=1
		return wheel	