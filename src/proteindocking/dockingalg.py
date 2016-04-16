#!/usr/bin/env python
# -*- coding: utf-8 -*-

from abc import ABCMeta, abstractmethod
from threading import Thread

class DockingAlgorithm(Thread):
    __metaclass__ = ABCMeta

    def __init__(self, docked_pair, stop_condition, selection, crossover):
        pass

    @abstractmethod
    def estimate_progress(self):
        return 0.0

    @abstractmethod
    def run(self):
        pass
