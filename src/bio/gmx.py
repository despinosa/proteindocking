import subprocess
import shutil
import os
import shlex
import sys
from tempfile import gettempdir

class gmx():

	TEMPDIR = gettempdir()
	ROOT = 'proteindocking'
	FILES = 'files'
	TMP = 'tmp'
	GMX_FILES = 'gmx_files'

	w = '10'
	h = '10'
	l = '10'

	protein_file = 'protein.pdb'
	#protein_hidro_file = 'protein_hidrogens.pdb'
	forcefield = 'charmm27'
	ligand_name = 'ligand'
	protein_ligand_box_file = 'conf_with_ligand.pdb'
	em_file = 'em.mdp'			
	topol_with_ligand_file = 'topol_with_ligand.top' 		

	@staticmethod
	def nspaces(n):
		cad = ''
		while(n):
			cad += ' '
			n -= 1
		return cad
	@staticmethod
	def add_molecule(line):
		pos = line.find('1')
		if pos == -1:
			raise Exception('Error')
		return (gmx.ligand_name+gmx.nspaces(pos-len(gmx.ligand_name)))+'1'
	@staticmethod	
	def ioFile():		
		with open(gmx.topol_with_ligand_file,'r') as in_file:
			buf = in_file.readlines()
		with open(gmx.topol_with_ligand_file,'w') as out_file:	
			allowed = 0		
			newfile = ''
			for i,line in enumerate(buf):			
				if '#include "'+gmx.forcefield+'.ff/forcefield.itp"' in line:
					line += '#include "'+gmx.ligand_name+'.itp"\n'				
				if '[ molecules ]' in line:
					allowed = 1
				if 'Protein' in line and allowed:					
					if i == len(buf) - 1:											
						line += gmx.add_molecule(line)
					elif '[' or ']' in line:						
						line += gmx.add_molecule(line)
				newfile += line
			out_file.write(newfile)		

	@staticmethod
	def process_topology():
		shutil.copyfile('topol.top',gmx.topol_with_ligand_file)
		gmx.ioFile()

	@staticmethod
	def generate_protein_topology():
		my_path = os.path.join(gmx.TEMPDIR, gmx.ROOT)
		if not os.path.exists(my_path):
			os.mkdir(my_path)
		my_path = os.path.join(my_path, gmx.FILES)
		if not os.path.exists(my_path):
			os.mkdir(my_path)
		os.chdir(my_path)
		cmd_protein_topology = "gmx pdb2gmx -ignh -f {} -ff {} -water none".format(gmx.protein_file,gmx.forcefield)
		try:
			p = subprocess.Popen(shlex.split(cmd_protein_topology), universal_newlines=True,
				                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out,err = p.communicate()
		except Exception:
			e = sys.exc_info()[1]
			print "Error: %s" % e
	@staticmethod
	def protein_ligand_box(): # thread_name, docked_pair): # no los usas
		########################################################################################################poner dockedpair
		cmd_protein_ligand_box = "gmx insert-molecules -f conf.gro -ci {}.pdb -o {} -box {} {} {} -nmol 1".format(gmx.ligand_name,gmx.protein_ligand_box_file,gmx.w,gmx.h,gmx.l)		
		try:
			p = subprocess.Popen(shlex.split(cmd_protein_ligand_box), universal_newlines=True,
				                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out,err = p.communicate()
		except Exception:
			e = sys.exc_info()[1]				
			print "Error: %s" % e

	@staticmethod
	def add_hydrogens():
		cmd_hydrogens = ["gmx grompp -f {} -o em_aux.tpr -c conf.gro".format(gmx.em_file),
					 "gmx trjconv -f conf.gro -o protein_hidrogens.pdb -s em_aux.tpr"] 
		try:
			n = len(cmd_hydrogens)
			i = 0
			while(i<n):
				p = subprocess.Popen(shlex.split(cmd_hydrogens[i]), universal_newlines=True)#,
					                     #stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				out,err = p.communicate()
				i += 1
		except Exception:
			e = sys.exc_info()[1]				
			print "Error: %s" % e


	@staticmethod
	def process_folders(thread_name):			
		os.chdir(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.TMP))

		if not os.path.exists(thread_name):			
			os.mkdir(thread_name)

		os.chdir(thread_name)

		if os.path.exists(gmx.GMX_FILES):
			shutil.rmtree(gmx.GMX_FILES)

		os.mkdir(gmx.GMX_FILES)
	
		shutil.copy(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.FILES, gmx.protein_file), gmx.GMX_FILES)
		#shutil.copy(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.FILES, gmx.protein_hidro_file), gmx.GMX_FILES)
		shutil.copy(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.FILES, '{}.itp'.format(gmx.ligand_name)), gmx.GMX_FILES)
		shutil.copy(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.FILES, '{}.pdb'.format(gmx.ligand_name)), gmx.GMX_FILES)
		shutil.copy(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.FILES, gmx.em_file), gmx.GMX_FILES)
		shutil.copy(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.FILES, gmx.topol_with_ligand_file), gmx.GMX_FILES)				
		shutil.copy(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.FILES, 'conf.gro'), gmx.GMX_FILES)				


		os.chdir(gmx.GMX_FILES)

	@staticmethod
	def calculate_fitness():
		energy = ''				
		cmd_fitness = ["gmx grompp -v -f {} -c {} -o em.tpr -p {}".format(gmx.em_file,gmx.protein_ligand_box_file,gmx.topol_with_ligand_file),
		   		   "gmx mdrun -v -s em.tpr"]
		try:
			n = len(cmd_fitness)
			i = 0
			while(i<n):
				p = subprocess.Popen(shlex.split(cmd_fitness[i]), universal_newlines=True,
			                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
		except Exception:
			e = sys.exc_info()[1]				
			print "Error: %s" % e
			print err
		return float(energy)

if __name__ == '__main__':
	#script_route = os.getcwd()
	gmx.generate_protein_topology()
	###############################gmx.add_hydrogens()  --opcional
	gmx.process_topology()
	gmx.process_folders('Hilo0')
	gmx.protein_ligand_box('Hilo0','../dockedpair')
	print gmx.calculate_fitness()


