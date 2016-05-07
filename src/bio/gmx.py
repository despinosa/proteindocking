import subprocess
import shutil
import os
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
    forcefield = 'charmm27'
    forcefield_atb = 'gromos54a7_atb.ff'
    topol_with_ligand_file = 'topol_with_ligand.top'   

    gmx_path = os.path.join(TEMPDIR, ROOT, TMP)        
    files_path = os.path.join(TEMPDIR, ROOT, FILES)      

    @staticmethod    
    def ioFile(ligand_name):        
        with open(gmx.topol_with_ligand_file,'r') as in_file:
            buf = in_file.readlines()
        with open(gmx.topol_with_ligand_file,'w') as out_file:    
            allowed = 0
            newfile = ''
            for i,line in enumerate(buf):            
                if '#include "'+gmx.forcefield+'.ff/forcefield.itp"' in line:
                    line += '#include "'+ligand_name+'.itp"\n'                
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
        cmd_center = ("gmx editconf -f {0} -c -o {1}".
                                format(molecule,molecule)
         try:
            p = subprocess.Popen(shlex.split(cmd_center), universal_newlines=True,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out,err = p.communicate()
            print out
            print err
        except Exception:
            e = sys.exc_info()[1]
            print "Error: %s" % e

    @staticmethod
    def process_topology(dp_object):
        os.chdir(gmx.files_path)
        shutil.copyfile('topol.top',gmx.topol_with_ligand_file)
        gmx.ioFile(dp_object.ligand_name)

    @staticmethod
    def generate_protein_topology(dp_object):                        
        os.chdir(gmx.files_path)
        cmd_protein_topology = ("gmx pdb2gmx -ignh -f {0} -ff {1} -water none".
                                format(dp_object.protein_file,gmx.forcefield))
        try:
            p = subprocess.Popen(shlex.split(cmd_protein_topology), universal_newlines=True,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out,err = p.communicate()
            print out
            print err
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
            while(i<n):
                p = subprocess.Popen(shlex.split(cmd_hydrogens[i]), universal_newlines=True,
                                         stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)                                                
                if i:
                    args = '0'                
                out,err = p.communicate(args)              
                print out
                print err                      
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
        energy = ''
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
            #obtiene la energia
            pos = err.find('Epot=')
            if pos != -1:
                i = pos + 5
                while(1):
                    if err[i] == 'F':
                        break
                    if err[i] != ' ':
                        energy += err[i]
                    i += 1
            ###################
            final_energy = float(energy)
        except Exception:
            e = sys.exc_info()[1]
            print "Error: %s" % e
            print err
            print thread_name
            p.kill()
            p.terminate()
            
        except ValueError:
            print energy
        return final_energy

