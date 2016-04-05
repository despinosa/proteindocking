from .molecule import Molecule

class Protein(Molecule):
    """docstring for Protein"""
    def __init__(self, pdb_id, pdb_file=None, cavities_pdb_file=None):
        super(Protein, self).__init__()
        self.pdb_id = pdb_id
        self.pdb_file = pdb_file
        self.cavities_pdb_file = cavities_pdb_file
        raise NotImplementedError
