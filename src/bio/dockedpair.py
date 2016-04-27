#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.PDB.DSSP import DSSP
from Bio.PDB.Model import Model
from Bio.PDB.PDBIO import Select, PDBIO as PDBOut # IO? solo escribe
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Vector import rotaxis2m, rotmat, Vector
from collections import namedtuple
from math import pi
from os import path, remove
import scipy as sp


class Model0Select(Select):
    def accept_model(self, model):
        return model.id == 0
        

class DockedPair(object):
    """docstring for DockedPair"""
    model0select = Model0Select()

    def __init__(self, main, arr):
        super(DockedPair, self).__init__()
        self.structure = main.original.copy()
        self.ligand = self.structure[0]['S']
        self.protein = self.structure[0]['P']
        self.cavities = self.structure[1]['C']
        self.decode(arr)

    def decode(self, arr):
        in_place = sp.array((0, 0, 0), 'f')
        rotation = rotaxis2m(arr[1], Vector(0, 0, 1))
        self.ligand.transform(rotation, in_place)
        rotation = rotaxis2m(arr[2], Vector(0, 1, 0))
        origin = self.cavities[int(round(arr[0]))]['R']
        self.ligand.transform(rotation, origin.coord)

    def free_energy(self, layer):
        out = PDBOut()
        out.set_structure(self.structure)
        path_ = 'tmp/{}/dockedpair.pdb'.format(layer.name)
        if path.exists(path_):
            out.save(path_, self.model0select)
        # pipes locos de gmx
        # parsear out de gmx y retorna energia (float)
        # borrar .pdb
        from random import random
        return random()

