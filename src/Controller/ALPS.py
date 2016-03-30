from Model.ALPSLayer import ALPSLayer
import threading

class ALPS:	
	AGE_GAP = 3
	def __init__(self,n_layers,n_individuals):
		self.layers = [ALPSLayer] * n_layers
		self.n_layers = n_layers
		self.n_individuals = n_individuals	

	@staticmethod
	def converge():
		return 0

	def generate_layers(self):
		max_age = 2
		for i in xrange(0,self.n_layers):			
			max_age = (max_age * 2) - 1
			self.layers[i] = ALPSLayer(max_age,self.n_individuals)
			self.layers[i].generate_population()	

	def runALPS(self):
		gen = 0
		while not ALPS.converge():
			for i in xrange(hola.n_layers-1,-1,-1):	
				upperLayer = None	
				if i < hola.n_layers-1:
					upperLayer = hola.layers[i+1]		
				if (i == 0) and (gen%hola.AGE_GAP):
					hola.layers[i].generate_population()	
					lowerLayer = None
				else:
					lowerLayer = hola.layers[i-1]		
					hola.layers[i].evolve(lowerLayer,upperLayer)
				print hola.layers[i]
				gen += 1		

hola = ALPS(12,10)
hola.generate_layers()
hola.runALPS()
