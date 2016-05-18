#!/usr/bin/env python
# -*- coding: utf-8 -*-

from alps.alps import ALPS
from bio.dockingproblem import DockingProblem
from alps.definitions.crossover import single_point
from bio.dockedpair import DockedPair
from Bio.PDB.PDBIO import Select
import Queue


class ALPSDocking(DockingProblem, ALPS):
    def __init__(self, ligand_path, protein_path, cavities_path, itp_path, preloaded_files_path,
                 forcefield, pop_size, mutate_rate, mating_rate, tourn_size,
                 stop_condition, elitism, aging_scheme, crossover=single_point,
                 n_parents=2, max_generations=1000, n_layers=10):
        super(ALPSDocking, self).__init__()
        DockingProblem.setup(self, ligand_path, protein_path, cavities_path,
                             itp_path, forcefield, preloaded_files_path)
        ALPS.setup(self, pop_size, mutate_rate, mating_rate, tourn_size,
                   stop_condition, elitism, aging_scheme, crossover, n_parents,
                   max_generations, n_layers)
        self.ex_queue = Queue.Queue()        

    def estimate_progress(self):
        return float(self.generation) / self.max_generations

    def check_errors(self):        
        for layer in self.layers:
            if not layer.ex_queue.empty():                
                return layer.ex_queue.get()

    def solve(self):
        self.start_layers() 
        self.join_layers()

    run = solve
    
    def _run_stdout(docking, output_path,pb_queue):    
        from datetime import datetime
        from sys import stdout
        from time import sleep    
        from os import path    

        docking.start()
        start = datetime.now()
        while docking.estimate_progress() < 1 - 1e-15:
            progress = docking.estimate_progress() * 100
            if pb_queue is not None:
                pb_queue.put(progress)
            stdout.write('\rprogress: \t{0:04.2f} %'.
                format(progress))
            ex = docking.check_errors()
            if ex is not None:
                raise ex
            stdout.flush()
            sleep(2)
        stdout.write('\n\n')
        docking.join()
        out_file = open(path.join(output_path,'info_best'),'w+')
        out_file.write('tiempo:\t{0}\n'.format(datetime.now()-start))
        pair = DockedPair(docking, docking.best)
        best_path = path.join(output_path,
                               'best_{0}_{1}.pdb'.format(docking.protein_filename.
                                        split('.')[0], docking.ligand_id))
        pair.to_file(best_path,
                     Select())
        out_file.write('mejor:\t{0}\n\n'.format(docking.best))
        out_file.close()      
        return best_path       

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
            ex = docking.check_errors()
            if ex is not None:
                raise ex
        pbar.finish()
        docking.join()
        out_file = open(path.join(output_path,'info_best'),'w+')
        out_file.write('tiempo:\t{0}\n'.format(datetime.now()-start))
        pair = DockedPair(docking, docking.best)
        best_path = path.join(output_path,
                               'best_{0}_{1}.pdb'.format(docking.protein_filename.
                                        split('.')[0], docking.ligand_id))
        pair.to_file(best_path,
                     Select())
        out_file.write('mejor:\t{0}\n\n'.format(docking.best))
        out_file.close()      
        return best_path

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

if __name__ == "__main__":
    from sys import argv
    from alps.definitions.agingscheme import fibonacci
    from alps.definitions.selection import enhanced
    from alps.definitions.stopcondition import conv_test, gen_limit
    from bio.dockedpair import DockedPair    
    
    (ligand_path, protein_path, cavities_path, itp_path, forcefield,
        output_path) = argv[1:7]
    forcefield = int(forcefield)
    fibo3 = lambda: fibonacci(3)
    # limited_conv = lambda alps: True if gen_limit(alps) else conv_test(alps)
    docking = ALPSDocking(ligand_path, protein_path, cavities_path, itp_path, 
                      output_path,
                      forcefield, 35, 0.1, 0.8, 5, gen_limit, enhanced,
                      fibo3, max_generations=33, n_layers=5)
    try:
        #docking._run_stdout(output_path,None)
        #docking._run_silent(output_path)
        docking._run_pbar(output_path)        
    except Exception as e:
        print e
