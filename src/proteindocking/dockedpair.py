#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.PDB.Structure import Structure

class DockedPair(Structure):
    """docstring for DockedPair"""
    def __init__(self, arg):
        super(DockedPair, self).__init__()
        self.arg = arg
