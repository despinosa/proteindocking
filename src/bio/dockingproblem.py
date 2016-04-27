#!/usr/bin/env python
# -*- coding: utf-8 -*-

from abc import ABCMeta, abstractmethod
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Model import Model
from Bio.PDB.PDBIO import Select, PDBIO as PDBOut # IO? solo escribe
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure
from collections import namedtuple
from dockedpair import DockedPair
from gmx import gmx
from math import pi
from threading import Thread
import scipy as sp


SecStructure = namedtuple('SecStructure', ('type_', 'first'))

class NonHOHSelect(Select):
    def accept_residue(self, residue):
        return residue.resname != 'HOH'

class DockingProblem(Thread):
    __metaclass__ = ABCMeta

    def setup(self, ligand_path, protein_path, cavities_path):
        def load_ligand(parser):
            ligand_model = parser.get_structure('ligand', ligand_path)[0]
            if len(ligand_model) != 1:
                raise ValueError(('el ligando ({}) no puede tener '
                                  'más de una cadena').format(ligand_path))
            ligand_chain = ligand_model.child_list[0]
            ligand_chain.id = 'S'
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
            for residue in cavities_chain.copy():
                cavities_chain.detach_child(residue.id)
                residue.id = (' ', residue.id[1], ' ')
                cavities_chain.add(residue)
            cavities_chain.id = 'C'
            return cavities_chain

        def encode(n_cavities):
            eps = 10**-15
            lower = sp.array((eps-0.5, 0, 0))
            upper = sp.array((n_cavities-eps-0.5, 2*pi, 2*pi))
            return lower, upper

        gmx.add_hydrogens()
        gmx.generate_protein_topology()
        gmx.process_topology()
        self.original = Structure('dockedpair')
        parser = PDBParser(PERMISSIVE=1)
        self.original.add(Model(0))
        self.original[0].add(load_protein(parser))
        self.original[0].add(load_ligand(parser))
        self.original.add(Model(1))
        cavities = load_cavities(parser)
        self.original[1].add(cavities)
        self.lower, self.upper = encode(len(cavities))

    def fitness(self, arr, layer):
        pair = DockedPair(self, arr)
        return pair.free_energy(layer)

    @abstractmethod
    def estimate_progress(self):
        return 0.0

    @abstractmethod
    def solve(self):
        pass

    run = solve


if __name__ == '__main__':
    from sys import argv
    from dockedpair import DockedPair
    problem = DockingProblem()
    problem.setup(*argv[1:4])
    arr = sp.array([7, 2.35, 1.88], 'f')
    problem.fitness(arr)
