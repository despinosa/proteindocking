#!/usr/bin/env python
# -*- coding: utf-8 -*-

from alps.alps import ALPS
from bio.dockingproblem import DockingProblem
from alps.definitions.agingscheme import linear

class ALPSDocking(DockingProblem, ALPS):
    def __init__(self, ligand_path, protein_path, cavities_path, itp_path,
                 forcefield, pop_size, mutate_rate, mating_rate, tourn_size,
                 stop_condition, elitism, crossover, n_parents=2,
                 max_generations=1000, aging_scheme=linear, n_layers=10):
        super(ALPSDocking, self).__init__()
        DockingProblem.setup(self, ligand_path, protein_path, cavities_path,
                             itp_path, forcefield)
        ALPS.setup(self, pop_size, mutate_rate, mating_rate, tourn_size,
                   stop_condition, elitism, crossover, n_parents,
                   max_generations, aging_scheme, n_layers)    
    
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
    from alps.definitions.stopcondition import conv_test, gen_limit
    from bio.dockedpair import DockedPair
    from Bio.PDB.PDBIO import Select

    ligand_path, protein_path, cavities_path, itp_path, forcefield = argv[1:6]
    limited_conv = lambda alps: True if gen_limit(alps) else conv_test(alps)
    docking = ALPSDocking(ligand_path, protein_path, cavities_path, itp_path,
                          forcefield, 15, 0.1, 0.8, 5, gen_limit, enhanced,
                          single_point, max_generations=111, n_layers=10)

    def run_pbar():
        from time import sleep
        from progressbar import ProgressBar, ReverseBar, ETA, Bar, Percentage    
        widgets = [Bar('>'), Percentage(),' ', ETA(), ' ', ReverseBar('<')]
        pbar = ProgressBar(widgets=widgets).start()  
        docking.start()
        start = datetime.now()  
        while docking.estimate_progress() < 1 - 1e-6:
            pbar.update(docking.estimate_progress() * 100)
        pbar.finish()
        docking.join()
        print 'tiempo:\t{0}\n'.format(datetime.now()-start)
        pair = DockedPair(docking, docking.best)
        pair.to_file(argv[6], Select())
        print 'mejor:\t{0}\n\n'.format(docking.best)

    def run_stdout():
        from sys import stdout
        from time import sleep
        docking.start()
        start = datetime.now()  
        while docking.estimate_progress() < 1 - 1e-6:
        stdout.write('\rprogreso:\t{0:04.2f} %'.
                format(docking.estimate_progress() * 100))
        stdout.flush()
        sleep(2)
        stdout.write('\n\n')
        pbar.finish()
        docking.join()
        print 'tiempo:\t{0}\n'.format(datetime.now()-start)
        pair = DockedPair(docking, docking.best)
        pair.to_file(argv[6], Select())
        print 'mejor:\t{0}\n\n'.format(docking.best)

    def run_silent():
        docking.start()
        docking.join()
        pair = DockedPair(docking, docking.best)
        pair.to_file(argv[6], Select())

