#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.PDB.DSSP import DSSP
from Bio.PDB.Model import Model
from Bio.PDB.PDBIO import Select, PDBIO as PDBOut # IO? solo escribe
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Vector import rotaxis2m, rotmat, Vector
from collections import namedtuple
from math import pi, cos, sin
from bisect import bisect
from os import mkdir, path, remove,rename
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
    all_ = Select()

    def __init__(self, main, arr):
        super(DockedPair, self).__init__()
        self.structure = main.original.copy()
        self.ligand_chain = self.structure[0]['Z']
        self.cavities_chain = self.structure[1]['C']
        self.main = main
        self.hash = arr.hash
        self.decode(arr)

    def decode(self, arr):
        cavity = self.cavities_chain[bisect(self.main.lise_rltt, arr[0])]['R']
        shift = cavity.occupancy * arr[3]
        origin = (np.array((shift * cos(arr[4]) * sin(arr[5]),
                            shift * sin(arr[4]) * sin(arr[5]),
                            shift * cos(arr[5])), 'f')
                  + cavity.coord)
        self.sqr_distance = sum(origin * origin)
        rotation = rotaxis2m(arr[1], Vector(0, 0, 1))
        self.ligand_chain.transform(rotation, origin)
        origin = np.array((0, 0, 0), 'f')
        rotation = rotaxis2m(arr[2], Vector(0, 1, 0))
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

