import subprocess
import shutil
import os
import re
import signal
import shlex
import sys
from tempfile import gettempdir
from threading import current_thread
from uuid import uuid4
from datetime import datetime

class gmx():

    TEMPDIR = gettempdir()        
    ROOT = 'proteindocking_{0}'.format(datetime.now().strftime('%Y%m%d%H%M%S%f'))
    FILES = 'files'
    TMP = 'tmp'
    GMX_FILES = 'gmx_files'    
    em_file = 'em.mdp'                
    CHARMM27, GROMOS54A7 = range(2)
    forcefields = {0:'charmm27',1:'gromos54a7_atb'}
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
                if forcefield_no == gmx.CHARMM27:        
                    if '#include "{0}.ff/forcefield.itp"'.format(gmx.forcefields[forcefield_no]) in line:
                        line += '#include "{0}.itp"\n'.format(ligand_name)                    
                elif forcefield_no == gmx.GROMOS54A7:                         
                    if '#include "./{0}.ff/forcefield.itp"'.format(gmx.forcefields[forcefield_no]) in line:                                                                                                                         
                        line += '#include "./{0}.ff/spc.itp"\n'.format(gmx.forcefields[forcefield_no])
                        line += '#include "./{0}.ff/ions.itp"\n'.format(gmx.forcefields[forcefield_no])
                        line += '#include "./{0}.itp"\n'.format(ligand_name)                        
                if '[ molecules ]' in line:
                    allowed = 1
                if 'protein' in line.lower() and allowed:
                    if i == len(buf) - 1:
                        line += ligand_name +' 1'
                    elif '[' or ']' or '\n' in line[i+1]:                        
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
            raise OSError('Checar GMX: No se genero el archivo de topologia de la proteina')
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
    def calculate_fitness(generation,hash_):        
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
            new_name = '{0}_{1}_{2}_{3}.pdb'.format(dockedpair.split('.')[0],generation,final_energy,hash_)
            if not os.path.isfile(new_name):
                os.rename(dockedpair,new_name)
        except ValueError:
            print str_energy.group()        
        except Exception:
            e = sys.exc_info()[1]    
            print "Error: %s" % e
            print out
            print err
            print thread_name                    
        return final_energy
