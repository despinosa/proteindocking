import subprocess
import shutil
import os
import re
import signal
import shlex
import sys
from tempfile import gettempdir
from threading import current_thread

class gmx():

    TEMPDIR = gettempdir()
    ROOT = 'proteindocking'
    FILES = 'files'
    TMP = 'tmp'
    GMX_FILES = 'gmx_files'    

    em_file = 'em.mdp'            
    forcefields = {1:'charmm27',2:'gromos54a7_atb'}
    topol_with_ligand_file = 'topol_with_ligand.top'   

    gmx_path = os.path.join(TEMPDIR, ROOT, TMP)        
    files_path = os.path.join(TEMPDIR, ROOT, FILES)      

    regexp_energy = re.compile('Epot=[ ]?[ -]?[\d]+[.]?[\d]+[eE]?[+-]?[\d]+')

    @staticmethod    
    def ioFile(ligand_name,forcefield_no):        
        with open(gmx.topol_with_ligand_file,'r') as in_file:
            buf = in_file.readlines()
        with open(gmx.topol_with_ligand_file,'w') as out_file:    
            allowed = 0
            newfile = ''
            for i,line in enumerate(buf):                    
                if forcefield_no == 1:        
                    if '#include "'+gmx.forcefields[forcefield_no]+'.ff/forcefield.itp"' in line:
                        line += '#include "'+ligand_name+'.itp"\n'
                elif forcefield_no == 2:
                    if '#include "'+gmx.forcefields[forcefield_no]+'.ff/forcefield.itp"' in line:     
                        line += ('#include "'+gmx.forcefields[forcefield_no]+'.ff/spc.itp"\n'
                        '#include "'+gmx.forcefields[forcefield_no]+'.ff/ions.itp"\n'
                        '#include "'+ligand_name+'.itp"\n')         
                if '[ molecules ]' in line:
                    allowed = 1
                if 'Protein' in line and allowed:
                    if i == len(buf) - 1:
                        line += ligand_name +' 1'
                    elif '[' or ']' in line:                        
                        line += ligand_name +' 1'
                newfile += line
            out_file.write(newfile)       
    @staticmethod
    def center_mol(molecule):
        os.chdir(gmx.files_path)
        molecule = molecule.split('.')[0]
        cmd_center = ("gmx editconf -f {0}.pdb -c -o {0}.pdb".format(molecule))
        try:
            p = subprocess.Popen(shlex.split(cmd_center), universal_newlines=True,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out,err = p.communicate()   
        except Exception:
            e = sys.exc_info()[1]
            print "Error: %s" % e

    @staticmethod
    def process_topology(dp_object):
        os.chdir(gmx.files_path)
        if not os.path.isfile('topol.top'):
            raise OSError('Checar GMX')
        shutil.copyfile('topol.top',gmx.topol_with_ligand_file)
        gmx.ioFile(dp_object.ligand_name,dp_object.forcefield)

    @staticmethod
    def generate_protein_topology(dp_object):                        
        os.chdir(gmx.files_path)
        cmd_protein_topology = ("gmx pdb2gmx -ignh -f {0} -ff {1} -water none -missing".
                                format(dp_object.protein_file,gmx.forcefields[dp_object.forcefield]))
        try:
            p = subprocess.Popen(shlex.split(cmd_protein_topology), universal_newlines=True,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out,err = p.communicate()            
        except Exception:
            e = sys.exc_info()[1]
            print "Error: %s" % e

    @staticmethod
    def add_hydrogens(dp_object):        
        os.chdir(gmx.files_path)
        cmd_hydrogens = [("gmx grompp -f {0} -o em_aux.tpr -c conf.gro".
                          format(gmx.em_file)),
                         ("gmx trjconv -f conf.gro -o {0} -s em_aux.tpr".
                          format(dp_object.protein_file))]
        try:
            n = len(cmd_hydrogens)
            i = 0
            args = None
            if not (os.path.isfile('conf.gro') or os.path.isfile(dp_object.protein_file)):
                raise OSError('Checar GMX')
            while(i<n):
                if i:
                    args = '0'
                    if not os.path.isfile('em_aux.tpr'):
                        raise OSError('Checar GMX')
                p = subprocess.Popen(shlex.split(cmd_hydrogens[i]), universal_newlines=True,
                                         stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)                                                                
                out,err = p.communicate(args)                                      
                p.stdin.close()                        
                i += 1
        except Exception:
            e = sys.exc_info()[1]                
            print "Error: %s" % e    
    @staticmethod
    def delete_files_in_path(dirPath):
        fileList = os.listdir(dirPath)
        for fileName in fileList:
             os.remove(os.path.join(dirPath,fileName))
                    
    @staticmethod
    def process_folders(dp_object):    
        try:        
            thread_name = current_thread().name                                            
            
            shutil.copy(os.path.join(gmx.files_path, dp_object.protein_file), gmx.gmx_path)        
            shutil.copy(os.path.join(gmx.files_path, '{0}.itp'.format(dp_object.ligand_name)), gmx.gmx_path)
            shutil.copy(os.path.join(gmx.files_path, '{0}.pdb'.format(dp_object.ligand_name)), gmx.gmx_path)
            shutil.copy(os.path.join(gmx.files_path, gmx.em_file), gmx.gmx_path)
            shutil.copy(os.path.join(gmx.files_path, gmx.topol_with_ligand_file), gmx.gmx_path)                
            shutil.copy(os.path.join(gmx.files_path, 'conf.gro'), gmx.gmx_path)                

        except Exception:
            e = sys.exc_info()[1]                
            print "Error: %s" % e
            print thread_name

    @staticmethod
    def calculate_fitness():        
        thread_name = current_thread().name    

        dockedpair = 'dockedpair_{0}.pdb'.format(thread_name)
        final_energy = float('inf')                
        cmd_fitness = [("gmx grompp -v -f {0} -c {1} -o em_{2}.tpr -p {3}".
                        format(gmx.em_file, dockedpair, thread_name,
                               gmx.topol_with_ligand_file)),
                       "gmx mdrun -v -s em_{0}.tpr".format(thread_name)]
        try:
            n = len(cmd_fitness)
            i = 0
            while(i<n):
                p = subprocess.Popen(shlex.split(cmd_fitness[i]),
                                     universal_newlines=True,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
                out,err = p.communicate()                              
                i += 1                            
                if p.returncode:                    
                    raise Exception(p.returncode)               
            str_energy = gmx.regexp_energy.search(err)                                  
            if str_energy:                                
                final_energy = float(str_energy.group().split('=')[-1].strip())                                
        except ValueError:
            print final_energy
        except Exception:
            e = sys.exc_info()[1]    
            print "Error: %s" % e
            print out
            print err
            print thread_name                    
        return final_energy
