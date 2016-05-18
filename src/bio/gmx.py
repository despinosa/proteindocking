import subprocess
import shutil
import re
import signal
import shlex
import sys
from os import path, environ, rename
from tempfile import gettempdir
from threading import current_thread
from uuid import uuid4
from datetime import datetime

class gmx():

    #Directorios
    TEMPDIR = gettempdir()        
    ROOT = 'proteindocking_{0}'.format(datetime.now().strftime('%Y%m%d%H%M%S%f'))
    FILES = 'files'
    TMP = 'tmp'
    GMX_FILES = 'gmx_files'                  
    gmx_path = path.join(TEMPDIR, ROOT, TMP)        
    files_path = path.join(TEMPDIR, ROOT, FILES)

    #Enum campos de fuerza
    CHARMM27, GROMOS54A7 = range(2)
    forcefields = {0:'charmm27',1:'gromos54a7_atb_pd'}

    #Archivos constantes para gmx
    topol_with_ligand_file = 'topol_with_ligand.top'   
    em_file = 'em.mdp'      
    topol_file = 'topol.top'
    confgro_file = 'conf.gro'
    em_aux_tpr_file = 'em_aux.tpr'         

    regexp_energy = re.compile('Epot=[ ]?[ -]?[\d]+[.]?[\d]+[eE]?[+-]?[\d]+')

    @staticmethod    
    def ioFile(topol_with_ligand_file_path,ligand_id,forcefield_no):        
        with open(topol_with_ligand_file_path,'r') as in_file:
            buf = in_file.readlines()
        with open(topol_with_ligand_file_path,'w') as out_file:    
            allowed = 0              
            newfile = ''
            for i,line in enumerate(buf):                                    
                if forcefield_no == gmx.CHARMM27:        
                    if '#include "{0}.ff/forcefield.itp"'.format(gmx.forcefields[forcefield_no]) in line:
                        ligand_file = path.join(gmx.files_path,ligand_id)
                        line += '#include "{0}.itp"\n'.format(ligand_file)                    
                elif forcefield_no == gmx.GROMOS54A7:                         
                    if '#include "{0}.ff/forcefield.itp"'.format(gmx.forcefields[forcefield_no]) in line:                                                                                                                         
                        line += '#include "{0}.ff/spc.itp"\n'.format(gmx.forcefields[forcefield_no])
                        line += '#include "{0}.ff/ions.itp"\n'.format(gmx.forcefields[forcefield_no])
                        line += '#include "{0}.itp"\n'.format(ligand_id)                        
                if '[ molecules ]' in line:
                    allowed = 1
                if 'protein' in line.lower() and allowed:
                    if i == len(buf) - 1:
                        line += ligand_id +' 1'
                    elif '[' or ']' or '\n' in line[i+1]:                        
                        line += ligand_id +' 1'                                
                newfile += line                
            out_file.write(newfile)     

    @staticmethod
    def center_mol(molecule):
        molecule = molecule.split('.')[0]+'.pdb'
        if not path.isfile(path.join(gmx.files_path,molecule)):
            raise OSError('Falta el archivo del ligando')        
        molecule_file = path.join(gmx.files_path,molecule).encode('string-escape')        
        cmd_center = ("gmx editconf -f {0} -c -o {1}".format(molecule_file,molecule_file.encode('unicode-escape')))
        
        p = subprocess.Popen(shlex.split(cmd_center), universal_newlines=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = p.communicate()           
        if p.returncode:                    
            raise Exception(str(p.returncode) + err)

    @staticmethod
    def process_topology(dp_object):        
        topol_file = path.join(gmx.files_path,gmx.topol_file)
        topol_with_ligand_file_path = path.join(gmx.files_path,gmx.topol_with_ligand_file) 
        if not path.isfile(topol_file):
            raise OSError('No se genero el archivo de topologia de la proteina')
        shutil.copyfile(topol_file,topol_with_ligand_file_path)
        gmx.ioFile(topol_with_ligand_file_path,dp_object.ligand_id,dp_object.forcefield)

    # Se usa string-escape para archivos de entrada *******************************************
    # unicode-escape para archivos de salida **************************************************
    @staticmethod
    def generate_protein_topology(dp_object):                          
        if not path.isfile(path.join(gmx.files_path,dp_object.protein_filename)):
            raise OSError('Falta el archivo de la proteina')

        protein_file = path.join(gmx.files_path,dp_object.protein_filename).encode('string-escape')              
        topol_out = path.join(gmx.files_path,gmx.topol_file).encode('unicode-escape')
        conf_out = path.join(gmx.files_path,gmx.confgro_file).encode('unicode-escape')
        posre_itp = path.join(gmx.files_path,'posre.itp').encode('unicode-escape')
        
        cmd_protein_topology = ("gmx pdb2gmx -ignh -f {0} -ff {1} -p {2} -o {3} -i {4} -water none -missing".
                                format(protein_file,gmx.forcefields[dp_object.forcefield],topol_out,conf_out,posre_itp))        
        p = subprocess.Popen(shlex.split(cmd_protein_topology), universal_newlines=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = p.communicate()                           
        if p.returncode:                    
            raise Exception(str(p.returncode) + err)       

    @staticmethod
    def add_hydrogens(dp_object):
        if not (path.isfile(path.join(gmx.files_path,gmx.confgro_file)) and path.isfile(path.join(gmx.files_path,gmx.em_file)) and path.isfile(path.join(gmx.files_path,gmx.topol_file))):
            raise OSError('Falta el archivo conf.gro o em.mdp')

        em_file_path = path.join(gmx.files_path,gmx.em_file).encode('string-escape')                     
        em_aux_out = path.join(gmx.files_path,gmx.em_aux_tpr_file).encode('unicode-escape')
        conf_file = path.join(gmx.files_path,gmx.confgro_file).encode('string-escape') 
        protein_file_out = path.join(gmx.files_path,dp_object.protein_filename).encode('unicode-escape')
        topol_file = path.join(gmx.files_path,gmx.topol_file).encode('string-escape')
        mdout_out = path.join(gmx.files_path,'mdout.mdp').encode('unicode-escape')
        cmd_hydrogens = [("gmx grompp -f {0} -o {1} -c {2} -p {3} -po {4}".
                          format(em_file_path,em_aux_out,conf_file,topol_file,mdout_out)),
                         ("gmx trjconv -f {0} -o {1} -s {2}".
                          format(conf_file,protein_file_out,em_aux_out.encode('string-escape')))]
        n = len(cmd_hydrogens)
        i = 0
        args = None        
        while(i<n):
            if i:
                args = '0'
                if not path.isfile(path.join(gmx.files_path,gmx.em_aux_tpr_file)):
                    raise OSError('Falta el archivo em_aux')
            p = subprocess.Popen(shlex.split(cmd_hydrogens[i]), universal_newlines=True,
                                     stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)                                                                
            out,err = p.communicate(args)                                
            if p.returncode:                    
                raise Exception(str(p.returncode) + err)
            p.stdin.close()                        
            i += 1                

    @staticmethod
    def preprocess(dp_object):                      
        gmx.center_mol(dp_object.ligand_id)
        gmx.generate_protein_topology(dp_object)        
        gmx.add_hydrogens(dp_object)
        gmx.generate_protein_topology(dp_object)
        gmx.process_topology(dp_object)
        environ['GMX_MAXBACKUP'] = '-1'

    @staticmethod
    def calculate_fitness(generation,hash_):        
        thread_name = current_thread().name 
        dockedpair = 'dockedpair_{0}.pdb'.format(thread_name)
        if not(path.isfile(path.join(gmx.files_path,gmx.em_file)) and path.isfile(path.join(gmx.files_path,gmx.topol_with_ligand_file)) and path.isfile(path.join(gmx.gmx_path,dockedpair))):
            raise OSError('Faltan archivos para el calculo de energia')
           
        em_file_path = path.join(gmx.files_path,gmx.em_file).encode('string-escape')
        topol_with_ligand_file_path = path.join(gmx.files_path,gmx.topol_with_ligand_file).encode('string-escape')
        em_thread_tpr_out = path.join(gmx.gmx_path,'em_{0}.tpr'.format(thread_name)).encode('unicode-escape')        
        dockedpair_file = path.join(gmx.gmx_path,dockedpair).encode('string-escape')
        traj_out = path.join(gmx.gmx_path,'traj.trr').encode('unicode-escape')
        mdout_out = path.join(gmx.gmx_path,'mdout.mdp').encode('unicode-escape')
        confoutgro_out = path.join(gmx.gmx_path,'confout.gro').encode('unicode-escape')
        ener_out = path.join(gmx.gmx_path,'ener.edr').encode('unicode-escape')
        mdlog_out = path.join(gmx.gmx_path,'md.log').encode('unicode-escape')
        
        final_energy = float('inf')                
        cmd_fitness = [("gmx grompp -v -f {0} -c {1} -o {2} -p {3} -po {4}".format(em_file_path, dockedpair_file, em_thread_tpr_out,topol_with_ligand_file_path,mdout_out)),
                       "gmx mdrun -v -s {0} -o {1} -c {2} -e {3} -g {4}".format(em_thread_tpr_out.encode('string-escape'), traj_out, confoutgro_out, ener_out, mdlog_out)]
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
                    raise Exception(str(p.returncode) + err)  
            str_energy = gmx.regexp_energy.search(err)
            if str_energy:
                final_energy = float(str_energy.group().split('=')[-1].strip())    
            new_name = '{0}_{1}_{2}_{3}.pdb'.format(dockedpair.split('.')[0],generation,final_energy,hash_)
            if not path.isfile(path.join(gmx.gmx_path,new_name)):
                rename(path.join(gmx.gmx_path,dockedpair),path.join(gmx.gmx_path,new_name))
        except ValueError:
            raise ValueError('Energia erronea')                            
        return final_energy
