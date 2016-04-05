#!/usr/bin/env python
# -*- coding: utf-8 -*-

from abc import ABCMeta, abstractmethod
from threading import Thread

class DockingAlg(object): #! Thread):
    __metaclass__ = ABCMeta
    @abstractmethod
    def estimate_progress(self):
        return 0.0
    @abstractmethod
    def run(self):
        pass
