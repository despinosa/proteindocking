from Model.ALPSLayer import ALPSLayer

class ALPS:	
	AGE_GAP = 3
	def __init__(self,n_layers,n_individuals):
		self.layers = [ALPSLayer] * n_layers
		self.n_layers = n_layers
		self.n_individuals = n_individuals

	def generate_layers(self):
		max_age = 2
		for i in range(0,self.n_layers):			
			max_age = (max_age * 2) - 1
			self.layers[i] = ALPSLayer(max_age,self.n_individuals)
			self.layers[i].generate_population()

hola = ALPS(12,10)
hola.generate_layers()
gen = 0
for i in range(hola.n_layers-1,-1,-1):		
	if (i == 0) and (gen%hola.AGE_GAP):
		hola.layers[i].generate_population()	
		lowerLayer = None
	else:
		lowerLayer = hola.layers[i-1]		
	hola.layers[i].evolve(lowerLayer)
gen += 1