#!/usr/bin/env python
# -*- coding: utf-8 -*-

from alps.alps import ALPS
from bio.dockingproblem import DockingProblem
from alps.definitions.crossover import single_point


class ALPSDocking(DockingProblem, ALPS):
    def __init__(self, ligand_path, protein_path, cavities_path, itp_path,
                 forcefield, pop_size, mutate_rate, mating_rate, tourn_size,
                 stop_condition, elitism, aging_scheme, crossover=single_point,
                 n_parents=2, max_generations=1000, n_layers=10):
        super(ALPSDocking, self).__init__()
        DockingProblem.setup(self, ligand_path, protein_path, cavities_path,
                             itp_path, forcefield)
        ALPS.setup(self, pop_size, mutate_rate, mating_rate, tourn_size,
                   stop_condition, elitism, aging_scheme, crossover, n_parents,
                   max_generations, n_layers)

    def estimate_progress(self):
        return float(self.generation) / self.max_generations

    def solve(self):
        self.start_layers()
        self.join_layers()

    run = solve


def _run_pbar(docking, output_path):
    from datetime import datetime
    from progressbar import ProgressBar, ReverseBar, ETA, Bar, Percentage     
    from os import path   
    widgets = [Bar('>'), Percentage(),' ', ETA(), ' ', ReverseBar('<')]
    pbar = ProgressBar(widgets=widgets).start()  
    docking.start()
    start = datetime.now()  
    while docking.estimate_progress() < 1 - 1e-15:
        pbar.update(docking.estimate_progress() * 100)
    pbar.finish()
    docking.join()
    out_file = open(path.join(output_path,'info_best'),'w+')
    out_file.write('tiempo:\t{0}\n'.format(datetime.now()-start))
    pair = DockedPair(docking, docking.best)
    pair.to_file(path.join(output_path,
                           'best_{0}_{1}.pdb'.format(docking.protein_filename.
                                    split('.')[0], docking.ligand_id)),
                 Select())
    out_file.write('mejor:\t{0}\n\n'.format(docking.best))
    out_file.close()

def _run_stdout(docking, output_path):
    from datetime import datetime
    from sys import stdout
    from time import sleep    
    from os import path
    docking.start()
    start = datetime.now()  
    while docking.estimate_progress() < 1 - 1e-15:
        stdout.write('\rprogreso:\t{0:04.2f} %'.
            format(docking.estimate_progress() * 100))
        stdout.flush()
        sleep(2)
    stdout.write('\n\n')
    docking.join()
    out_file = open(path.join(output_path,'info_best'),'w+')
    out_file.write('tiempo:\t{0}\n'.format(datetime.now()-start))
    pair = DockedPair(docking, docking.best)
    pair.to_file(path.join(output_path,
                           'best_{0}_{1}.pdb'.format(docking.protein_filename.
                                    split('.')[0], docking.ligand_id)),
                 Select())
    out_file.write('mejor:\t{0}\n\n'.format(docking.best))
    out_file.close()

def _run_silent(docking, output_path):
    from datetime import datetime    
    from os import path
    docking.start()
    start = datetime.now()  
    docking.join()
    out_file = open(path.join(output_path,'info_best'),'w+')
    out_file.write('tiempo:\t{0}\n'.format(datetime.now()-start))
    pair = DockedPair(docking, docking.best)
    pair.to_file(path.join(output_path,
                           'best_{0}_{1}.pdb'.format(docking.protein_filename.
                                    split('.')[0], docking.ligand_id)),
                 Select())
    out_file.write('mejor:\t{0}\n\n'.format(docking.best))
    out_file.close()

if __name__ == '__main__':
    from sys import argv
    from alps.definitions.agingscheme import fibonacci
    from alps.definitions.selection import enhanced
    from alps.definitions.stopcondition import conv_test, gen_limit
    from bio.dockedpair import DockedPair
    from Bio.PDB.PDBIO import Select    

    (ligand_path, protein_path, cavities_path, itp_path, forcefield,
            output_path) = argv[1:7]
    forcefield = int(forcefield)
    fibo3 = lambda: fibonacci(3)
    # limited_conv = lambda alps: True if gen_limit(alps) else conv_test(alps)
    docking = ALPSDocking(ligand_path, protein_path, cavities_path, itp_path,
                          forcefield, 35, 0.1, 0.8, 5, gen_limit, enhanced,
                          fibo3, max_generations=48, n_layers=5)
    _run_pbar(docking, output_path)
