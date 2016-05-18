from sys import stdout
from os import path    
from datetime import datetime
from time import sleep
from alpsdocking import ALPSDocking
from alps.alps import ALPS
from bio.dockingproblem import DockingProblem
from alps.definitions.crossover import single_point
from alps.definitions.agingscheme import fibonacci
from alps.definitions.selection import enhanced
from alps.definitions.stopcondition import gen_limit  
from pymol import cmd
from threading import Thread

class ALPSMain(Thread):
    def __init__(self,queue,protein_path,ligand_path,itp_path,cavities_path,output_path,forcefield,files_path):
        super(ALPSMain, self).__init__()
        self.queue = queue
        self.protein_path = str(protein_path)
        self.ligand_path = str(ligand_path)
        self.itp_path = str(itp_path)
        self.cavities_path = str(cavities_path)
        self.output_path = output_path
        self.forcefield = forcefield        
        self.FILES_PATH = files_path    

    def run(self):        
        fibo3 = lambda: fibonacci(3)    
        docking = ALPSDocking(self.ligand_path, self.protein_path, self.cavities_path, self.itp_path, self.FILES_PATH,
                              self.forcefield, 10, 0.1, 0.8, 5, gen_limit, enhanced,
                              fibo3, max_generations=5, n_layers=5)
        pair_file = docking._run_stdout(self.output_path,self.queue)               
        self.queue.put("Docking complete.\nPDB file: {0}".format(pair_file))
        cmd.load(pair_file)