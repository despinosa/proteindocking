from Model.Chromosome import Chromosome

class ALPSLayer:	
	def __init__(self,max_age,n_individuals):
		self.max_age = max_age
		self.n_individuals = n_individuals		
		self.population = list();
	@classmethod
	def constructor_knownBest(self,max_age,n_individuals,best):
		self.max_age = max_age
		self.n_individuals = n_individuals
		self.best = best
		self.population = list()
	
	def generate_population(self):
		for i in range(0,hola.n_individuals):					
			hola.population.append(Chromosome())	
	    
hola = ALPSLayer(10,10)
hola.generate_population();
hola.population.sort()
for i in range(0,hola.n_individuals):	
	print hola.population[i],
	print hola.population[i].score