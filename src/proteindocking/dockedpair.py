#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.PDB.DSSP import DSSP
from Bio.PDB.PDBIO import Select, PDBIO as PDBOut # IO? solo escribe
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure


class DockedPair(Structure):
    """docstring for DockedPair"""

    class NonHOHSelect(Select):
        def accept_residue(self, residue):
            return residue.get_name() != 'HOH'

    class Model0Select(Select):
        def accept_model(self, model):
            return model.get_id() == 0


    def __init__(self, ligand_path, protein_path, cavities_path):
        def load_ligand(parser):
            ligand_model = parser.get_structure('ligand', ligand_path)[0]
            if len(ligand_model) != 1:
                raise ValueError(('el ligando ({}) no puede tener '
                                  'más de una cadena').format(ligand_path))
            dssp = DSSP(ligand_model, ligand_path)
            ligand_chain = ligand_model.child_list[0]
            for res_info in dssp:
                ligand_chain[res_info[0]].xtra['sec_struct'] = res_info[2]
            ligand_chain.id = 'C'
            return ligand_chain

        def load_protein(parser):
            protein_struct = parser.get_structure('protein', protein_path)
            out = PDBOut()
            out.set_structure(protein_struct)
            out.save(protein_path, NonHOHSelect())
            protein_struct = parser.get_structure('protein', protein_path)
            protein_model = protein_struct[0]
            if len(protein_model) != 1:
                raise ValueError(('la proteína ({}) no puede tener '
                                  'más de una cadena').format(protein_path))
            protein_chain = protein_model.child_list[0]
            protein_chain.id = 'C'
            return protein_chain

        def load_cavities(parser):
            cavities_model = parser.get_structure('cavities', cavities_path)[0]
            if len(cavities_model) != 1:
                raise ValueError(('las cavidades ({}) deben estar '
                                  'en una sola cadena').format(cavities_path))
            cavities_chain = cavities_model.child_list[0]
            cavities_chain.id = 'C'
            return cavities_chain

        super(DockedPair, self).__init__(self, 'dockedpair')
        parser = PDBParser(PERMISSIVE=1)
        self.ligand = load_ligand(parser)
        self.add(Model(0))
        self[0].add(self.ligand)
        self.protein = load_protein(parser)
        self[0].add(self.protein)
        self.cavities = load_cavities(parser)
        self.add(Model(1))
        self[1].add(self.cavities)


if __name__ == '__main__':
    from sys import argv
    pair = DockedPair(*argv[1:4])
