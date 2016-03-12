from Functions.Functions import Functions
import random

class Chromosome:
	L_INDIVIDUAL = 10
	def __init__(self):
		self.score = random.random();
		self.age = 1;
		self.chromosome = Functions.random_chromosome(self.L_INDIVIDUAL)

	@classmethod
	def constructor_knownAge(self,age):
		self.score = 0
		self.age = age			
	
	def __str__(self):
		s = ""
		for i in range(0,self.L_INDIVIDUAL):
			s += str(self.chromosome[i])
		return s

	def __cmp__(self,other):		
		if self.score < other.score:
			return -1
		elif self.score > other.score:
			return 1
		else:
			return 0
	def incrementAge():
		self.chromosome.age += 1
	# def mutate():

	# def scoring():

	
