#!/usr/bin/env python
# -*- coding: utf-8 -*-

from alpslayer import ALPSLayer
from dockingalg import DockingAlg
from ..util import Singleton

class ALPS(DockingAlg):
    __metaclass__ = Singleton

    def __init__(self, problem, stop_condition, selection, crossover):
        pass

    def estimate_progress(self):
        pass
