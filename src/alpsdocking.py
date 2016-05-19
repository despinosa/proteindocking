#!/usr/bin/env python
# -*- coding: utf-8 -*-

from alps.alps import ALPS
from alps.definitions.crossover import single_point
from alps.definitions.agingscheme import fibonacci
from alps.definitions.selection import enhanced
from alps.definitions.stopcondition import gen_limit
from bio.dockedpair import DockedPair
from bio.dockingproblem import DockingProblem
from Bio.PDB.PDBIO import Select
from os import path
import Queue

class ALPSDocking(DockingProblem, ALPS):
    def __init__(self, ligand_path, protein_path, cavities_path, itp_path, preloaded_files_path,
                 forcefield):        
        def config_file():
            f = open(path.join(preloaded_files_path,'files','config'), 'r')
            self.config_args = []
            for i,line in enumerate(f):
                self.config_args.append(line.split('=')[-1]) 
            if(len(self.config_args) != 7):
                raise Exception('Archivo de configuracion erroneo.')            
        super(ALPSDocking, self).__init__()        
        self.ex_queue = Queue.Queue() 
        config_file()
        aging_scheme_factor = int(self.config_args[0])
        pop_size = int(self.config_args[1])
        mutate_rate = float(self.config_args[2])
        mating_rate = float(self.config_args[3])
        tourn_size = int(self.config_args[4])
        max_generations = int(self.config_args[5])
        n_layers = int(self.config_args[6])
        fibo = lambda: fibonacci(aging_scheme_factor)

        DockingProblem.setup(self, ligand_path, protein_path, cavities_path,
                             itp_path, forcefield, preloaded_files_path)
        ALPS.setup(self, pop_size, mutate_rate, mating_rate, tourn_size,
                   gen_limit, enhanced, fibo, single_point, max_generations, n_layers)
               

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
    from sys import argv, exc_info
    from alps.definitions.agingscheme import fibonacci
    from alps.definitions.selection import enhanced
    from alps.definitions.stopcondition import gen_limit
    from bio.dockedpair import DockedPair
    from datetime import datetime
    from traceback import format_exc    
    import logging        
    
    (ligand_path, protein_path, cavities_path, itp_path, forcefield,
        output_path) = argv[1:7]
    logger = None
    try:
        docking = ALPSDocking(ligand_path, protein_path, cavities_path, itp_path, output_path, int(forcefield))
        #docking._run_stdout(output_path,None)
        #docking._run_silent(output_path)
        docking._run_pbar(output_path)        
    except Exception as e:
        if logger == None:
            logging.basicConfig(filename='exceptions_{0}.log'.format(datetime.now().strftime('%Y%m%d%H%M%S%f')),level=logging.DEBUG)
            logger = logging.getLogger('Protein docking')
        logger.info(str(e)+'\n'+format_exc()+'\n'+str(exc_info()[0]))
