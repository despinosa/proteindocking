import collections
from SecondaryStructure import SecondaryStructure

class Molecule:
	def __init__(self,pdb_file_location,structure):
		self.pdb_file_location = pdb_file_location
		self.structure = structure
