from sys import argv, stdout
from os import chdir, path    
from datetime import datetime
from time import sleep
from alpsdocking import ALPSDocking
from alps.alps import ALPS
from bio.dockingproblem import DockingProblem
from alps.definitions.crossover import single_point
from alps.definitions.agingscheme import fibonacci
from alps.definitions.selection import enhanced
from alps.definitions.stopcondition import gen_limit
from bio.dockedpair import DockedPair
from Bio.PDB.PDBIO import Select    
from pymol import cmd
from threading import Thread

class ALPSMain(Thread):
    def __init__(self,queue,protein_path,ligand_path,itp_path,cavities_path,output_path,forcefield,files_path):
        super(ALPSMain, self).__init__()
        self.queue = queue
        self.protein_path = protein_path
        self.ligand_path = ligand_path
        self.itp_path = itp_path
        self.cavities_path = cavities_path
        self.output_path = output_path
        self.forcefield = forcefield        
        self.FILES_PATH = files_path

    def dock(self,docking):                            
        docking.start()
        start = datetime.now()              
        while docking.estimate_progress() < 1 - 1e-15:                                      
            stdout.write('\rprogress: \t{0:04.2f} %'.
                format(docking.estimate_progress()*100))
            stdout.flush()
            sleep(2)                        
        stdout.write('\n\n')
        docking.join()
        out_file = open(path.join(self.output_path,'info_best'),'w+')
        out_file.write('tiempo:\t{0}\n'.format(datetime.now()-start))
        pair = DockedPair(docking, docking.best)
        self.pair_file = path.join(self.output_path,'best_{0}_{1}.pdb'.format(docking.protein_filename.split('.')[0], docking.ligand_id))
        pair.to_file(self.pair_file,Select())
        out_file.write('mejor:\t{0}\n\n'.format(docking.best))
        out_file.close()        
        cmd.load(self.pair_file)

    def run(self):        
        fibo3 = lambda: fibonacci(3)    
        docking = ALPSDocking(self.ligand_path, self.protein_path, self.cavities_path, self.itp_path, self.FILES_PATH,
                              self.forcefield, 35, 0.1, 0.8, 5, gen_limit, enhanced,
                              fibo3, max_generations=5, n_layers=5)
        self.dock(docking)
        self.queue.put("Docking complete.\nPDB file: {0}".format(self.pair_file)) #escribe en la cola para que el proceso de progress bar termine