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
from os import chdir,mkdir,path,remove,getcwd,makedirs,environ
from re import match, search
from shutil import rmtree, copy
from threading import Thread
import scipy as sp


SecStructure = namedtuple('SecStructure', ('type_', 'first'))

class NonHOHSelect(Select):
    def accept_residue(self, residue):
        return residue.resname != 'HOH'

class DockingProblem(Thread):
    __metaclass__ = ABCMeta

    def setup(self, ligand_path, protein_path, cavities_path,itp_path,forcefield):
        def load_ligand(parser):
            ligand_model = parser.get_structure('ligand',
                                                self.new_ligand_path)[0]
            if len(ligand_model) != 1:
                raise ValueError(('el ligando ({0}) no puede tener más de una '
                                  'cadena').format(self.new_ligand_path))
            ligand_chain = ligand_model.child_list[0]
            ligand_chain.id = 'S'
            return ligand_chain

        def load_protein(parser):
            protein_struct = parser.get_structure('protein',
                                                  self.new_protein_path)
            out = PDBOut()
            out.set_structure(protein_struct)
            out.save(self.new_protein_path, NonHOHSelect())
            protein_struct = parser.get_structure('protein',
                                                  self.new_protein_path)
            protein_model = protein_struct[0]
            if len(protein_model) != 1:
                raise ValueError(('la proteína ({0}) no puede tener más de una'
                                  ' cadena').format(self.new_protein_path))
            protein_chain = protein_model.child_list[0]
            protein_chain.id = 'P'
            return protein_chain

        def load_cavities(parser):
            cavities_model = parser.get_structure('cavities',
                                                  self.new_cavities_path)[0]
            if len(cavities_model) != 1:
                raise ValueError(('las cavidades ({0}) deben estar en una sola'
                                  ' cadena').format(self.new_cavities_path))
            cavities_chain = cavities_model.child_list[0]
            for residue in cavities_chain.copy():
                cavities_chain.detach_child(residue.id)
                residue.id = (' ', residue.id[1], ' ')
                cavities_chain.add(residue)
            cavities_chain.id = 'C'
            return cavities_chain

        def load_folders():            
            if path.exists(path.join(gmx.TEMPDIR,gmx.ROOT)):
                rmtree(path.join(gmx.TEMPDIR,gmx.ROOT),ignore_errors=True)
            mkdir(path.join(gmx.TEMPDIR,gmx.ROOT))
            mkdir(gmx.files_path)
            mkdir(gmx.gmx_path)

        def load_files():                                    
            copy(ligand_path, path.join(gmx.files_path,
                                        self.ligand_name+'.pdb'))
            copy(itp_path, path.join(gmx.files_path,
                                        self.ligand_name+'.itp'))
            copy(protein_path, gmx.files_path)
            copy(cavities_path, gmx.files_path)
            copy(path.join('files', gmx.em_file),
                 gmx.files_path)

        def id_from_itp(itp_file):
            moltype = r'\s*\[\s*(moleculetype|moltype)\s*\]'
            comment = r'\s*;'
            id_ = r'\s*(\w+)'
            for line in itp_file:
                if match(moltype, line.lower()): break
            else:
                raise ValueError('no se encuentra el identificador en {}'.
                                        format(itp_file.name))
            line = next(itp_file)
            while match(comment, line): line = next(itp_file)
            try:
                return search(id_, line).group(1)
            except AttributeError:
                raise ValueError('no se encuentra el identificador en {}'.
                                        format(itp_file.name))

        load_folders()
        self.protein_file = path.split(protein_path)[1]
        self.ligand_name = id_from_itp(open(itp_path, 'r'))
        self.cavities_file = path.split(cavities_path)[1]
        self.forcefield = int(forcefield)
        load_files()
        self.new_protein_path = path.join(gmx.files_path,
                                          self.protein_file)
        self.new_ligand_path = path.join(gmx.files_path,
                                         '{0}.pdb'.format(self.ligand_name))
        self.new_cavities_path = path.join(gmx.files_path,
                                           self.cavities_file)
        gmx.center_mol(self.protein_file)
        gmx.center_mol(self.ligand_name)
        gmx.generate_protein_topology(self)               
        remove(self.new_protein_path)
        gmx.add_hydrogens(self)
        gmx.generate_protein_topology(self)
        gmx.process_topology(self)
        gmx.process_folders(self)        
        environ['GMX_MAXBACKUP'] = '-1'
        chdir(gmx.gmx_path)

        self.original = Structure('dockedpair')
        parser = PDBParser(PERMISSIVE=1)
        self.original.add(Model(0))
        self.protein = load_protein(parser)
        self.original[0].add(self.protein)
        self.ligand = load_ligand(parser)
        self.original[0].add(self.ligand)
        self.original.add(Model(1))
        self.cavities = load_cavities(parser)
        self.original[1].add(self.cavities)
        self.encode()

    def encode(self):
        self.lise_rltt = map(lambda cav: cav['R'].bfactor, self.cavities)
        for i in xrange(1, len(self.lise_rltt)):
            self.lise_rltt[i] += self.lise_rltt[i-1]
        lise_max = self.lise_rltt.pop()
        self.lower = sp.array((     0.0,  0.0,  0.0, 0.0, 2*pi, 2*pi), 'f')
        self.upper = sp.array((lise_max, 2*pi, 2*pi, 1.0, 2*pi, 2*pi), 'f')

    def fitness(self, arr):
        pair = DockedPair(self, arr)
        return pair.free_energy()

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
