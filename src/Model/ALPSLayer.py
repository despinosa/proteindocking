from Model.Chromosome import Chromosome
from Functions.Functions import Functions
import random

CROSSOVER_POINT = Chromosome.L_CHROMOSOME/2	
N_ELITISM = 2 #define cuantos individuos pasaran directamente a la siguiente generacion

class ALPSLayer:	
	def __init__(self,max_age,n_individuals):
		self.max_age = max_age
		self.n_individuals = n_individuals		
		self.population = [Chromosome] * n_individuals				
	
	def generate_population(self):
		for i in range(0,self.n_individuals):					
			self.population[i] = Chromosome()	

	def elitism(self,nextGeneration):
		for i in range(0,N_ELITISM):
			self.population[i].incrementAge()
			nextGeneration[i] = self.population[i]

	def crossover(self,parent1,parent2,previous_layer): # cruzamiento
		global CROSSOVER_POINT		
		if parent1 < self.n_individuals:
			newChromosome1=self.population[parent1].chromosome[0:CROSSOVER_POINT]	
		elif parent1 >= self.n_individuals:			
			newChromosome1=previous_layer.population[parent1%previous_layer.n_individuals].chromosome[0:CROSSOVER_POINT]

		if parent2 < self.n_individuals:
			newChromosome1.extend(self.population[parent2].chromosome[CROSSOVER_POINT:Chromosome.L_CHROMOSOME])
			newChromosome2=self.population[parent2].chromosome[0:CROSSOVER_POINT]
		elif parent2 >= self.n_individuals:
			parent2 %= previous_layer.n_individuals
			newChromosome1.extend(previous_layer.population[parent2].chromosome[CROSSOVER_POINT:Chromosome.L_CHROMOSOME])		
			newChromosome2=previous_layer.population[parent2].chromosome[0:CROSSOVER_POINT]
			
		if parent1 < self.n_individuals:
			newChromosome2.extend(self.population[parent1].chromosome[CROSSOVER_POINT:Chromosome.L_CHROMOSOME])	
		elif parent1 >= self.n_individuals:
			parent1 %= previous_layer.n_individuals
			newChromosome2.extend(previous_layer.population[parent1].chromosome[CROSSOVER_POINT:Chromosome.L_CHROMOSOME])	
					
		newAge = self.population[parent1].age if self.population[parent1].age >= self.population[parent2].age else self.population[parent2].age

		descendant1 = Chromosome.constructor_crossover(newAge,newChromosome1)
		descendant2 = Chromosome.constructor_crossover(newAge,newChromosome2)
		return descendant1,descendant2

	def evolve(self,previous_layer):
		self.population.sort()		
		nextGeneration = self.population[:]
		self.elitism(nextGeneration)
		roulette = Functions.create_wheel(self,previous_layer)					
		for i in range(0,(self.n_individuals-N_ELITISM)/2):
			parent1 = random.choice(roulette)
			parent2 = random.choice(roulette)					
			descendant1, descendant2 = self.crossover(parent1,parent2,previous_layer)	
			descendant1.mutate()		
			descendant2.mutate()
			nextGeneration[N_ELITISM + 2*i] = descendant1
			nextGeneration[(N_ELITISM + 1) + 2*i] = descendant2
		self.population = nextGeneration[:]
