from Model.ALPSLayer import ALPSLayer
import threading

class ALPS:		
	def __init__(self,n_layers,n_individuals,n_elitism,n_parents,tourn_size,maximum_generations,age_gap):
		self.layers = [ALPSLayer] * n_layers
		self.n_layers = n_layers
		self.n_individuals = n_individuals
		self.n_parents = n_parents
		self.n_elitism = n_elitism
		self.tourn_size = tourn_size		
		self.maximum_generations = maximum_generations
		self.gen = 0
		self.age_gap = age_gap

	def converge(self):
		convergence = 0
		if convergence or (self.gen >= self.maximum_generations):
			return 1
		return 0

	def generate_layers(self):
		max_age = 2
		for i in xrange(0,self.n_layers):
			if i == self.n_layers-1:
				max_age = float('inf')
			else:
				max_age = (max_age * 2) - 1 #AGE GAP				
			self.layers[i] = ALPSLayer(max_age,self.n_individuals,self.n_elitism,self.n_parents,self.tourn_size)			
		self.layers[0].generate_population()	

	def assign_layers(self):
		for i, lay in enumerate(self.layers[1:]):
			lay.previous_layer = self.layers[i-1]
		for i, lay in enumerate(self.layers[:-1]):
			lay.next_layer = self.layers[i+1]

	def run(self):		
		while not self.converge():
			for i in xrange(self.n_layers-1,-1,-1):							
				if (i == 0) and (self.gen%self.age_gap == 0):
					self.layers[i].generate_population()										
				self.layers[i].evolve()		
				#print self.gen						
				self.gen += 1
				
				
	def printpruebas(self):
		for i in xrange(0,self.n_layers):
			if self.layers[i].populated:
				print "LAYER" + str(i)
				print "MAX_AGE" + str(self.layers[i].max_age)
				print "INDIVIDUOS" + str(len(self.layers[i].population))
 				print self.layers[i]


hola = ALPS(10,200,5,2,5,100,3)
hola.generate_layers()
hola.assign_layers()
hola.run()
#hola.printpruebas()
