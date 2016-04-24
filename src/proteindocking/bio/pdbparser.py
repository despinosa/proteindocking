from Tkinter import * 
from tkFileDialog import *
import os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.DSSP import DSSP

class PDBParser(PDBParser):
    """docstring for PDBParser"""
    OPTIONS = {}
    OPTIONS['defaultextension'] = '.pdb'
    OPTIONS['filetypes'] = [('PDB files','.pdb'),('ENT files','.ent'),('All files', '.*')]

    #def __init__(self, file_loc, type_):
    def __init__(self):
        super(PDBParser,self).__init__()
        #self.file_loc = file_loc
        #elf.type_ = type_

        
    def load(self):
        os.chdir('../../..')                
        os.environ["PATH"] += os.pathsep + os.getcwd()
        root = Tk()
        root.withdraw()
        file_loc = askopenfilename(**self.OPTIONS)                
        structure = self.get_structure("test",file_loc)        
        model = structure[0]
        dssp = DSSP(model, file_loc)                       
        a_key = list(dssp)                
        # for r in dssp:
        #     print(r)
        # print("Handled %i residues" % len(dssp))
        # print(sorted(dssp))
        # if ('A', 1) in dssp:
        #     print(dssp[('A', 1)])
        #     print(structure[0]['A'][1].xtra)
        # Secondary structure
       #print(''.join(item[1] for item in dssp))
        for item in a_key:            
            if item[2] == 'E' or item[2] == 'H':
                print item  
        pass

    def load_cavities():
        pass

    def load_nm():
        pass


p = PDBParser()
p.load()
