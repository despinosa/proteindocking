from abc import ABCMeta, abstractmethod
from threading import Thread

class DockingAlg(Thread):
    __metaclass__ = ABCMeta
    @abstractmethod
    def estimate_progress(self):
        pass
    @abstractmethod
    def run(self):
        pass
