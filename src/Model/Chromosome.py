from Functions.Functions import Functions
import random

class Chromosome:
	L_CHROMOSOME = 50
	MUTATION_PROBABILITY = 0.7

	def __init__(self):
		self.score = random.random();
		self.age = 1;
		self.chromosome = Functions.random_chromosome_numpy(self.L_CHROMOSOME)
		print self.chromosome

	@classmethod
	def constructor_crossover(self,age,chromosome):
		self = Chromosome()
		self.score = random.random();
		self.age = age			
		self.chromosome = chromosome
		return self
	
	def __str__(self):		
		return str(self.chromosome)

	def __cmp__(self,other):				
		if self.score < other.score:
			return -1
		elif self.score > other.score:
			return 1
		else:
			return 0	
	def mutate(self):
		if random.random() < self.MUTATION_PROBABILITY:
			self.chromosome[random.randint(0,self.L_CHROMOSOME-1)] = Functions.random_float_value()		

	def incrementAge(self):
		self.age += 1

	# def scoring():

	
