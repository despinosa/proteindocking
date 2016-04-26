import subprocess
import shutil
import os
import shlex
import sys

class gmx():

	protein_file = '1ubq.pdb'
	forcefield = 'charmm27'
	ligand_name = 'zinc_11909586'
	protein_ligand_box_file = 'conf_with_ligand10x10x10.pdb'
	em_file = 'em.mdp'		
	files_route = '../tmp/'
	docked_files_route = '../../tmp/'
	docked_pair = 'dockedpair.pdb'

	protein_topology = "gmx pdb2gmx -ignh -f "+ files_route + protein_file + " -ff " + forcefield + " -water none"

	

	cmd_hydrogens = ["gmx grompp -f em.mdp -o em_aux.tpr -c conf.gro",
					 "gmx trjconv -f conf.gro -o conf_final.pdb -s em_aux.tpr"] 

	cmd_fitness = ["gmx grompp -v -f em.mdp -c conf_with_ligand10x10x10.pdb -o em.tpr -p topol_with_ligand.top",
		   		   "gmx mdrun -v -s em.tpr"]	

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
	def generate_protein_topology():
		try:
			p = subprocess.Popen(shlex.split(gmx.protein_topology), universal_newlines=True,
				                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out,err = p.communicate()
		except Exception:
			e = sys.exc_info()[1]				
			print "Error: %s" % e
	@staticmethod
	def protein_ligand_box(thread_name):
		cmd_protein_ligand_box = "gmx insert-molecules -f conf.gro -ci "+gmx.docked_files_route+thread_name+gmx.docked_pair+" -o conf_with_ligand10x10x10.pdb -box 10 10 10 -nmol 1"
		try:
			p = subprocess.Popen(shlex.split(cmd_protein_ligand_box), universal_newlines=True,
				                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out,err = p.communicate()
		except Exception:
			e = sys.exc_info()[1]				
			print "Error: %s" % e

	@staticmethod
	def add_hydrogens():
		try:
			n = len(gmx.cmd_hydrogens)
			i = 0
			while(i<n):
				p = subprocess.Popen(shlex.split(gmx.cmd_hydrogens[i]), universal_newlines=True,
					                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				out,err = p.communicate()
				i += 1
		except Exception:
			e = sys.exc_info()[1]				
			print "Error: %s" % e

	@staticmethod	
	def ioFile():		
		with open('topol_with_ligand.top','r') as in_file:
			buf = in_file.readlines()
		with open('topol_with_ligand.top','w') as out_file:	
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
		shutil.copyfile('topol.top','topol_with_ligand.top')
		gmx.ioFile()	

	@staticmethod
	def process_folders(thread_name):
		tmpfolder = thread_name

		if os.path.exists(tmpfolder):
			shutil.rmtree(tmpfolder)

		os.mkdir(tmpfolder)
	
		shutil.copy(gmx.files_route+gmx.protein_file,tmpfolder)
		shutil.copy(gmx.files_route+gmx.ligand_name+'.itp',tmpfolder)
		shutil.copy(gmx.files_route+gmx.ligand_name+'.pdb',tmpfolder)
		shutil.copy(gmx.files_route+gmx.em_file,tmpfolder)
		shutil.copy(gmx.files_route+gmx.protein_ligand_box_file,tmpfolder)

	@staticmethod
	def calculate_fitness(thread_name):				
		energy = ''				
		tmpfolder = thread_name		
		try:
			n = len(gmx.cmd_fitness)
			i = 0
			while(i<n):								
				p = subprocess.Popen(shlex.split(gmx.cmd_fitness[i]), universal_newlines=True,
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
			os.chdir('..')			
			if os.path.exists(tmpfolder):
				shutil.rmtree(tmpfolder)			
		except Exception:
			e = sys.exc_info()[1]				
			print "Error: %s" % e
			print err

		return energy

#Ejemplo
os.chdir(gmx.files_route) # cambia al directorio donde estan los archivos
gmx.generate_protein_topology() #genera la topologia de la proteina
gmx.process_topology() # genera la topologia de la proteina con el ligando
os.chdir('../bio') # cambia de vuelta al directorio donde esta el codigo

#fitness
gmx.process_folders('HiloPrueba') 
os.chdir('HiloPrueba')
print gmx.calculate_fitness('HiloPrueba')