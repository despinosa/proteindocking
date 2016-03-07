from Molecule import Molecule

class Ligand(Molecule):
	def __init__(self,flexible_points):
		self.fps = flexible_points;