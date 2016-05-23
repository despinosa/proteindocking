#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.PDB.DSSP import DSSP
from Bio.PDB.Model import Model
from Bio.PDB.PDBIO import Select, PDBIO as PDBOut # IO? solo escribe
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Vector import rotaxis2m, rotmat, Vector
from collections import namedtuple
from math import pi, cos, sin, exp
from bisect import bisect
from os import mkdir, path, remove ,rename
from tempfile import gettempdir
from threading import current_thread
from gmx import gmx
import numpy as np


class Model0Select(Select):
    def accept_model(self, model):
        return model.id == 0


class DockedPair(object):
    """docstring for DockedPair"""
    model0 = Model0Select()

    def __init__(self, main, arr):
        super(DockedPair, self).__init__()
        self.structure = main.original.copy()
        self.ligand_chain = self.structure[0]['Z']
        self.main = main
        self.hash = arr.hash
        self.decode(arr)

    def decode(self, arr):
        rotation = rotaxis2m(arr[0], Vector(0, 0, 1))
        origin = np.array((0, 0, 0), 'f')
        self.ligand_chain.transform(rotation, origin)
        rotation = rotaxis2m(arr[1], Vector(0, 1, 0))
        self.shift = arr[2]
        origin = (np.array((self.shift * sin(arr[3]) * cos(arr[4]),
                            self.shift * sin(arr[3]) * sin(arr[4]),
                            self.shift * cos(arr[3])), 'f'))
        self.ligand_chain.transform(rotation, origin)

    def to_file(self, pdb_path, select=model0):
        out = PDBOut()
        out.set_structure(self.structure)        
        out.save(pdb_path, select) #, self.model0select)  

    def free_energy(self):
        pdb_path = path.join(gmx.gmx_path,
                             'dockedpair_{0}.pdb'.format(current_thread().name))
        self.to_file(pdb_path)
        return gmx.calculate_fitness(self.main.generation, self.hash)

