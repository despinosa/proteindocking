#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.PDB.DSSP import DSSP
from Bio.PDB.Model import Model
from Bio.PDB.PDBIO import Select, PDBIO as PDBOut # IO? solo escribe
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Vector import rotmat, Vector
from collections import namedtuple
from math import pi
from scipy import identity

class DockedPair(Structure):
    """docstring for DockedPair"""

    class NonHOHSelect(Select):
        def accept_residue(self, residue):
            return residue.get_name() != 'HOH'

    class Model0Select(Select):
        def accept_model(self, model):
            return model.get_id() == 0

    SecStructure = namedtuple('SecStructure', ('type_', 'first'))

    def __init__(self, ligand_path=None, protein_path=None, cavities_path=None):
        def load_ligand(parser):
            ligand_model = parser.get_structure('ligand', ligand_path)[0]
            if len(ligand_model) != 1:
                raise ValueError(('el ligando ({}) no puede tener '
                                  'más de una cadena').format(ligand_path))
            dssp = DSSP(ligand_model, ligand_path)
            ligand_chain = ligand_model.child_list[0]
            for res_info in dssp:
                ligand_chain[res_info[0]].xtra['sec_structure'] = res_info[2]
            ligand_chain.id = 'L'
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
            protein_chain.id = 'P'
            return protein_chain

        def load_cavities(parser):
            cavities_model = parser.get_structure('cavities', cavities_path)[0]
            if len(cavities_model) != 1:
                raise ValueError(('las cavidades ({}) deben estar '
                                  'en una sola cadena').format(cavities_path))
            cavities_chain = cavities_model.child_list[0]
            cavities_chain.id = 'C'
            return cavities_chain

        def encode(ligand_chain, n_cavities):
            eps = 10**-15
            encoding = {0 : (eps - 0.5, n_cavities - 0.5 - eps)}
            type_ = ligand_chain[0].xtra['sec_structure']
            sec_structures = []
            first = 0
            for i, residue in enumerate(ligand_chain):
                if residue.xtra['sec_structure'] != type_:
                    encoding[len(encoding)] = (0, 2*pi)
                    sec_structure = SecStructure(first, type_)
                    sec_structures.append(sec_structure)
                    first = i
                type_ = residue.xtra['sec_structure']
            return encoding, sec_structures

        super(DockedPair, self).__init__(self) # , 'dockedpair')
        if (ligand_path, protein_path, cavities_path) is not (None, None, None):
            parser = PDBParser(PERMISSIVE=1)
            self.add(Model(0))
            self[0].add(load_ligand(parser))
            self[0].add(load_protein(parser))
            self.add(Model(1))
            self.__class__.cavities = load_cavities(parser)
            self[1].add(self.cavities)
            self.__class__.encoding, self.__class__.sec_structures = encode(
                self[0]['L'], len(self.cavities)
            )
            self.__class__.state0 = self.copy()
        else:
            self = self.state0

    def decode(self, array):
        def deform(angles):
            origin = -self.ligand[self.sec_structures[0].first]['CA']
            for i, angle in enumerate(angles):
                rotation = rotmat(angle, origin.normalized())
                for j in xrange(self.sec_structures[i].first,
                                len(self.ligand)):
                    self.ligand[j].transform(rotation, origin)
                if i+1 < len(angles):
                    origin = (-self.ligand[self.sec_structures[i+1].first]
                              ['CA'].coord)

        def translate(origin):
            self.transform(identity(3), -origin)

        self = self.state0
        deform(array[1:])
        translate(self.cavities[int(round(array[0]))])


if __name__ == '__main__':
    from sys import argv
    pair = DockedPair(*argv[1:4])
