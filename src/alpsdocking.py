#!/usr/bin/env python
# -*- coding: utf-8 -*-

from alps.alps import ALPS
from bio.dockingproblem import DockingProblem
from os import mkdir, path
from tempfile import gettempdir

class ALPSDocking(DockingProblem, ALPS):
    def __init__(self, ligand_path, protein_path, cavities_path, # docking
                 pop_size, mutate_rate, mating_rate, tourn_size, # alps
                 stop_condition, elitism, crossover):
        super(ALPSDocking, self).__init__()
        ALPS.setup(self, pop_size, mutate_rate, mating_rate, tourn_size,
                   stop_condition, elitism, crossover)
        DockingProblem.setup(self, ligand_path, protein_path, cavities_path)
        root = path.join(gettempdir(), 'tmp')
        if not path.exists(root):
            mkdir(root)
        for layer in self.layers:
            new_folder = path.join(root, layer.name)
            if not path.exists(new_folder):
                mkdir(new_folder)


    def estimate_progress(self):
        return self.generation / self.max_generations

    def solve(self):
        self.start_layers()
        self.join_layers()

    run = solve


if __name__ == '__main__':
    from datetime import datetime
    from sys import argv
    from alps.definitions.crossover import single_point
    from alps.definitions.selection import enhanced
    from alps.definitions.stopcondition import gen_limit

    ligand_path, protein_path, cavities_path = argv[1:4]
    docking = ALPSDocking(ligand_path, protein_path, cavities_path,
                          10, 0.1, 0.8, 5, gen_limit, enhanced,
                          single_point)
    docking.start()
    start = datetime.now()
    docking.join()
    print 'time: {}\n'.format(datetime.now()-start)
    print 'best: {}\n\n'.format(docking.layers[-1].population[0])
