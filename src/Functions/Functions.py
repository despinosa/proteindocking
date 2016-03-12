import random



class Functions:
	@staticmethod
	def random_chromosome(L_INDIVIDUAL):
		chromosome=[]
		for i in range(0,L_INDIVIDUAL):
			if random.random()<0.7:
				chromosome.append(0)
			else:
				chromosome.append(1)
		return chromosome