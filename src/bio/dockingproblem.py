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
from os import chdir,mkdir,path,remove,getcwd,makedirs
import shutil


SecStructure = namedtuple('SecStructure', ('type_', 'first'))

class NonHOHSelect(Select):
    def accept_residue(self, residue):
        return residue.resname != 'HOH'

class DockingProblem(Thread):
    __metaclass__ = ABCMeta

    def setup(self, ligand_path, protein_path, cavities_path,itp_path):
        def load_ligand(parser):
            ligand_model = parser.get_structure('ligand', self.new_ligand_path)[0]
            if len(ligand_model) != 1:
                raise ValueError(('el ligando ({}) no puede tener '
                                  'más de una cadena').format(self.new_ligand_path))
            ligand_chain = ligand_model.child_list[0]
            ligand_chain.id = 'S'
            return ligand_chain

        def load_protein(parser):
            protein_struct = parser.get_structure('protein', self.new_protein_path)
            out = PDBOut()
            out.set_structure(protein_struct)
            out.save(self.new_protein_path, NonHOHSelect())
            protein_struct = parser.get_structure('protein', self.new_protein_path)
            protein_model = protein_struct[0]
            if len(protein_model) != 1:
                raise ValueError(('la proteína ({}) no puede tener '
                                  'más de una cadena').format(self.new_protein_path))
            protein_chain = protein_model.child_list[0]
            protein_chain.id = 'P'
            return protein_chain

        def load_cavities(parser):
            cavities_model = parser.get_structure('cavities', self.new_cavities_path)[0]
            if len(cavities_model) != 1:
                raise ValueError(('las cavidades ({}) deben estar '
                                  'en una sola cadena').format(self.new_cavities_path))
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

        def load_folders():
            currentwd = getcwd()           
            if path.exists(path.join(gmx.TEMPDIR,gmx.ROOT)):
                shutil.rmtree(path.join(gmx.TEMPDIR,gmx.ROOT),ignore_errors=True)            
            mkdir(path.join(gmx.TEMPDIR,gmx.ROOT))            
            mkdir(path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES))
            mkdir(path.join(gmx.TEMPDIR,gmx.ROOT,gmx.TMP))
            root = path.join(gmx.TEMPDIR,gmx.ROOT,gmx.TMP)
            for layer in self.layers:
                new_folder = path.join(root, layer.name,gmx.GMX_FILES)
                if not path.exists(new_folder):
                    makedirs(new_folder)
            chdir(currentwd)       

        def load_files():                                    
            shutil.copy(ligand_path,path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES))  
            shutil.copy(itp_path,path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES))
            shutil.copy(protein_path,path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES))              
            shutil.copy(cavities_path,path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES)) 
            shutil.copy('tmp/'+gmx.em_file,path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES)) 

        load_folders()
        self.protein_file = path.split(protein_path)[1]        
        self.ligand_name = path.split(ligand_path)[1].split('.')[0]        
        self.cavities_file = path.split(cavities_path)[1]
        load_files()
        self.new_protein_path = path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES,self.protein_file)
        self.new_ligand_path = path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES,'{}.pdb'.format(self.ligand_name))
        self.new_cavities_path = path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES,self.cavities_file)
        gmx.generate_protein_topology(self)   
        currentwd = getcwd()
        chdir(path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES))        
        remove(self.protein_file)
        gmx.add_hydrogens(self)
        gmx.generate_protein_topology(self)
        gmx.process_topology(self)
        chdir(currentwd)

        self.original = Structure('dockedpair')
        parser = PDBParser(PERMISSIVE=1)
        self.original.add(Model(0))
        self.original[0].add(load_protein(parser))
        self.original[0].add(load_ligand(parser))
        self.original.add(Model(1))
        cavities = load_cavities(parser)
        self.original[1].add(cavities)
        self.lower, self.upper = encode(len(cavities))

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
