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
from math import exp, pi
from os import mkdir,path,remove,getcwd,makedirs,environ
from re import match, search
from shutil import rmtree, copy, copytree
from threading import Thread
import numpy as np


SecStructure = namedtuple('SecStructure', ('type_', 'first'))

class NonHOHSelect(Select):
    def accept_residue(self, residue):
        return residue.resname != 'HOH'

class DockingProblem(Thread):
    __metaclass__ = ABCMeta

    def setup(self, ligand_path, protein_path, cavities_path, itp_path,
              forcefield,preloaded_files):
        def load_data():
            parser = PDBParser(PERMISSIVE=1)
            protein_struct = parser.get_structure('protein',
                                                  self.protein_path)
            out = PDBOut()
            out.set_structure(protein_struct)
            out.save(self.protein_path, NonHOHSelect())
            self.protein_model = parser.get_structure('protein',
                                                  self.protein_path)[0]
            ligand_model = parser.get_structure('ligand',
                                                self.ligand_path)[0]
            if len(ligand_model) != 1:
                raise ValueError(('el ligando ({0}) no puede tener más de una '
                                  'cadena').format(self.ligand_path))
            self.ligand_chain = ligand_model.child_list[0]
            self.ligand_chain.id = 'Z'
            cavities_model = parser.get_structure('cavities',
                                                  self.cavities_path)[0]
            if len(cavities_model) != 1:
                raise ValueError(('las cavidades ({0}) deben estar en una sola'
                                  ' cadena').format(self.cavities_path))
            self.cavities_chain = cavities_model.child_list[0]
            for residue in self.cavities_chain.copy():
                self.cavities_chain.detach_child(residue.id)
                residue.id = (' ', residue.id[1], ' ')
                self.cavities_chain.add(residue)
            self.cavities_chain.id = 'C'
            self.original = Structure('dockedpair')
            self.original.add(self.protein_model)
            self.original[0].add(self.ligand_chain)
            self.original.add(Model(1))
            self.original[1].add(self.cavities_chain)

        def prepare_wdtree():
            if path.exists(path.join(gmx.TEMPDIR,gmx.ROOT)):
                rmtree(path.join(gmx.TEMPDIR,gmx.ROOT), ignore_errors=True)
            mkdir(path.join(gmx.TEMPDIR,gmx.ROOT))
            mkdir(gmx.files_path)
            mkdir(gmx.gmx_path)

        def prepare_wdfiles():                                    
            copy(ligand_path, path.join(gmx.files_path,
                                        '{0}.pdb'.format(self.ligand_id)))
            copy(itp_path, path.join(gmx.files_path,
                                     '{0}.itp'.format(self.ligand_id)))
            copy(protein_path, gmx.files_path)
            copy(cavities_path, gmx.files_path)                        
            copy(path.join(preloaded_files,'files', gmx.em_file),
                 gmx.files_path)
            if(self.forcefield == gmx.GROMOS54A7):
                environ['GMXDATA'] = path.join(gmx.files_path)
                environ['GMXLIB'] = path.join(gmx.files_path)
                copytree(path.join(preloaded_files,'files', gmx.forcefields[gmx.GROMOS54A7]+'.ff'),
                         path.join(gmx.files_path,gmx.forcefields[gmx.GROMOS54A7]+'.ff'))                
                copy(path.join(preloaded_files,'files','residuetypes.dat'), gmx.files_path)
                copy(path.join(preloaded_files,'files','elements.dat'), gmx.files_path)
                copy(path.join(preloaded_files,'files','xlateat.dat'), gmx.files_path)
                copy(path.join(preloaded_files,'files','specbond.dat'), gmx.files_path)
            self.protein_path = path.join(gmx.files_path, self.protein_filename)
            self.ligand_path = path.join(gmx.files_path,
                                         '{0}.pdb'.format(self.ligand_id))
            self.cavities_path = path.join(gmx.files_path,
                                           self.cavities_filename)

        def id_from_itp(itp_file):
            moltype = r'\s*\[\s*(moleculetype|moltype)\s*\]'
            comment = r'\s*;'
            id_ = r'\s*(\w+)'
            not_found = ValueError('no se encuentra el identificador en {0}'.
                                       format(itp_file.name))
            for line in itp_file:
                if match(moltype, line.lower()): break
            else:
                raise not_found
            line = next(itp_file)
            while match(comment, line): line = next(itp_file)
            try:
                return search(id_, line).group(1)
            except AttributeError:
                raise not_found

        prepare_wdtree()
        self.protein_filename = path.split(protein_path)[-1]
        with open(itp_path, 'r') as itp_file:
            self.ligand_id = id_from_itp(itp_file)
        self.cavities_filename = path.split(cavities_path)[-1]
        self.forcefield = forcefield
        prepare_wdfiles()
        gmx.preprocess(self)                      
        load_data()
        self.encode()        

    def encode(self):
        """Codifica el problema en un arreglo de longitud 6.

        Se construye una _ruleta_ de selección para determinar la
        cavidad a probar en una solución propuesta. También se definen
        los límites superior e inferior de cada posición del arreglo.

        arr[0]  ~   Valor del tiro de la _ruleta_ de cavidades.
        arr[1:2]~   Ángulos de rotación del ligando.
        arr[3:5]~   Desfase del ligando en coordenadas esféricas respecto
                    del centro de la cavidad.

        """

        self.lise_rltt = map(lambda cav: exp(cav['R'].bfactor),
                             self.cavities_chain)
        for i in xrange(1, len(self.lise_rltt)):
            self.lise_rltt[i] += self.lise_rltt[i-1]
        lise_max = self.lise_rltt.pop()
        self.lower = np.array((     0.0,  0.0,  0.0, 0.0,  0.0,  0.0), 'f')
        self.upper = np.array((lise_max, 2*pi, 2*pi, 1.0, 2*pi, 2*pi), 'f')
        self.span = self.upper - self.lower

    def fitness(self, arr):        
        pair = DockedPair(self, arr)
        return pair.free_energy() + exp(10 * pair.shift / pair.cavity.bfactor)

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
    arr = np.array([7, 2.35, 1.88], 'f')
    problem.fitness(arr)
