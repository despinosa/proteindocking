#!/usr/bin/env python
# -*- coding: utf-8 -*-

from alps.alps import ALPS
from bio.dockingproblem import DockingProblem
from alps.definitions.agingscheme import linear

class ALPSDocking(DockingProblem, ALPS):
    def __init__(self, ligand_path, protein_path, cavities_path, itp_path, # docking
                 pop_size, mutate_rate, mating_rate, tourn_size, # alps
                 stop_condition, elitism, crossover, n_parents=2,
                 max_generations=50, aging_scheme=linear, n_layers=5):
        super(ALPSDocking, self).__init__()
        ALPS.setup(self, pop_size, mutate_rate, mating_rate, tourn_size,
                   stop_condition, elitism, crossover, n_parents,
                   max_generations, aging_scheme, n_layers)
        DockingProblem.setup(self, ligand_path, protein_path, cavities_path,
                             itp_path)        
    
    def estimate_progress(self):
        return float(self.generation) / self.max_generations

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
    from bio.dockedpair import DockedPair
    from Bio.PDB.PDBIO import Select
    from sys import stdout
    from time import sleep

    ligand_path, protein_path, cavities_path, itp_path = argv[1:5]
    docking = ALPSDocking(ligand_path, protein_path, cavities_path, itp_path,
                          10, 0.1, 0.8, 5, gen_limit, enhanced,
                          single_point)
    docking.start()
    start = datetime.now()
    while docking.estimate_progress() < 1 - 1e-6:
        stdout.write('\rprogress: {:04.2f} %'.format(docking.estimate_progress() * 100))
        stdout.flush()
        sleep(2)
    stdout.write('\n\n')
    docking.join()
    print 'time: {}\n'.format(datetime.now()-start)
    best = docking.layers[0].population[0]
    for i in xrange(1, max(1, len(docking.layers)-2)):
        try:
            best = min(best, docking.layers[i].population[0])
        except IndexError:
            pass
    pair = DockedPair(docking, best)
    pair.to_file(argv[5], Select())
    print 'best: {}\n\n'.format(best)
