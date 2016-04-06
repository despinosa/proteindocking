from Model.Chromosome import Chromosome
from Functions.Functions import Functions
import random
import threading
import numpy as np


class ALPSLayer(threading.Thread):	
	def __init__(self,max_age,n_individuals,n_elitism,n_parents,tourn_size):
		super(ALPSLayer, self).__init__()
		self.max_age = max_age
		self.n_individuals = n_individuals		
		self.population = list()	
		self.n_elitism = n_elitism
		self.n_parents = n_parents
		self.tourn_size = tourn_size
		self.nextGeneration = None
		self.populated = 0

	def __str__(self):
		s = ""
		for i in xrange(0,self.n_individuals):
			s += str(self.population[i]) + "\n"
			s += "AGE:" + str(self.population[i].age)
		return s
	
	def generate_population(self):		
		self.population = [Chromosome] * self.n_individuals	
		for i in xrange(0,self.n_individuals):					
			self.population[i] = Chromosome()	
		self.populated = 1

	def elitism(self):
		for i in xrange(0,self.n_elitism):
			self.population[i].incrementAge()
			self.nextGeneration[i] = self.population[i]

	def crossover_roulette(self,parent1,parent2,previous_layer,upper_layer): # cruzamiento		
		CROSSOVER_POINT = random.randint(1,Chromosome.L_CHROMOSOME - 1)		
		if parent1 < self.n_individuals: # se checa de que poblacion viene el parent escogido ya que la ruleta regresa un numero de 0 a 2*n_individuals - 1
			newChromosome1=self.population[parent1].chromosome[0:CROSSOVER_POINT]	
		elif parent1 >= self.n_individuals:			
			newChromosome1=previous_layer.population[parent1%previous_layer.n_individuals].chromosome[0:CROSSOVER_POINT]

		if parent2 < self.n_individuals:
			newChromosome1 = np.append(newChromosome1, self.population[parent2].chromosome[CROSSOVER_POINT:Chromosome.L_CHROMOSOME], 0)
			newChromosome2=self.population[parent2].chromosome[0:CROSSOVER_POINT]
		elif parent2 >= self.n_individuals:
			parent2 %= previous_layer.n_individuals
			newChromosome1 = np.append(newChromosome1, previous_layer.population[parent2].chromosome[CROSSOVER_POINT:Chromosome.L_CHROMOSOME], 0)		
			newChromosome2=previous_layer.population[parent2].chromosome[0:CROSSOVER_POINT]
			
		if parent1 < self.n_individuals:
			newChromosome2 = np.append(newChromosome2, self.population[parent1].chromosome[CROSSOVER_POINT:Chromosome.L_CHROMOSOME], 0)	
		elif parent1 >= self.n_individuals:
			parent1 %= previous_layer.n_individuals
			newChromosome2 = np.append(newChromosome2, previous_layer.population[parent1].chromosome[CROSSOVER_POINT:Chromosome.L_CHROMOSOME], 0)	
					
		newAge = self.population[parent1].age if self.population[parent1].age >= self.population[parent2].age else self.population[parent2].age

		descendant1 = Chromosome.constructor_crossover(newAge,newChromosome1)
		descendant2 = Chromosome.constructor_crossover(newAge,newChromosome2)									
		descendant1.incrementAge()
		descendant2.incrementAge()
		return descendant1,descendant2	

	def crossover(self,parents):		
		CROSSOVER_POINT = random.randint(1,Chromosome.L_CHROMOSOME - 1)
		newChromosome1 = parents[0].chromosome[0:CROSSOVER_POINT]
		newChromosome1 = np.append(newChromosome1,parents[1].chromosome[CROSSOVER_POINT:Chromosome.L_CHROMOSOME],0)
		newChromosome2 = parents[1].chromosome[0:CROSSOVER_POINT]
		newChromosome2 = np.append(newChromosome2,parents[0].chromosome[CROSSOVER_POINT:Chromosome.L_CHROMOSOME],0)		
		return newChromosome1,newChromosome2

	def generateOffspring(self,previous_layer,upper_layer):
		parents = self.tourney(previous_layer)	
		newAge = parents[0].age if parents[0].age >= parents[1].age else parents[1].age
		newChromosome1,newChromosome2 = self.crossover(parents)
		descendant1 = Chromosome.constructor_crossover(newAge,newChromosome1)
		descendant2 = Chromosome.constructor_crossover(newAge,newChromosome2)									
		descendant1.incrementAge()
		descendant2.incrementAge()
		if(descendant1.age >= self.max_age):
			upper_layer.population.append(descendant1)
			upper_layer.populated = 1			
		if(descendant2.age >= self.max_age):
			upper_layer.population.append(descendant2)	
			upper_layer.populated = 1			
		return descendant1,descendant2				


	def tourney(self, previous_layer):
		if previous_layer == None:
			source = self.population
		else:
			source = previous_layer.population+self.population
		source_size = len(source)
		tournament = random.sample(source, min(self.tourn_size,source_size))				
		tournament.sort()						
		return tournament[:self.n_parents]		

	def buildNextGeneration(self,i,descendant1,descendant2):
		self.nextGeneration[self.n_elitism + 2*i] = descendant1
		self.nextGeneration[(self.n_elitism + 1) + 2*i] = descendant2

	def replacement(self):
		self.population = self.nextGeneration[:]		

	def evolve(self,previous_layer,upper_layer):
		try:
			self.population.sort()		
			self.nextGeneration = self.population[:]
			self.elitism()
			# roulette = Functions.create_wheel(self,previous_layer)
			for i in xrange(0,(len(self.population)-self.n_elitism)/2):									
				descendant1, descendant2 = self.generateOffspring(previous_layer,upper_layer)
				descendant1.mutate()		
				descendant2.mutate()
				self.buildNextGeneration(i,descendant1,descendant2)
			################################################################################################################	
			self.replacement()
		except Exception as ex:
			template = "An exception of type {0} occured. Arguments:\n{1!r}"
			message = template.format(type(ex).__name__, ex.args)
			print message	
