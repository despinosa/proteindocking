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

	CMD = ["gmx pdb2gmx -ignh -f "+ protein_file + " -ff " + forcefield + " -water none",
		   "gmx grompp -v -f em.mdp -c conf_with_ligand10x10x10.pdb -o em.tpr -p topol_with_ligand.top",
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
		return (gmxfitness.ligand_name+gmxfitness.nspaces(pos-len(gmxfitness.ligand_name)))+'1'

	@staticmethod	
	def ioFile():		
		with open('topol_with_ligand.top','r') as in_file:
			buf = in_file.readlines()
		with open('topol_with_ligand.top','w') as out_file:	
			allowed = 0		
			newfile = ''
			for i,line in enumerate(buf):			
				if '#include "'+gmxfitness.forcefield+'.ff/forcefield.itp"' in line:
					line += '#include "'+gmxfitness.ligand_name+'.itp"\n'				
				if '[ molecules ]' in line:
					allowed = 1
				if 'Protein' in line and allowed:					
					if i == len(buf) - 1:											
						line += gmxfitness.add_molecule(line)
					elif '[' or ']' in line:						
						line += gmxfitness.add_molecule(line)
				newfile += line
			out_file.write(newfile)

	@staticmethod
	def process_topology():		
		shutil.copyfile('topol.top','topol_with_ligand.top')
		gmxfitness.ioFile()	

	
	@staticmethod
	def calculate_fitness(thread_name):				
		energy = ''		
		tmpfolder = thread_name

		if os.path.exists(tmpfolder):
			shutil.rmtree(tmpfolder)

		os.mkdir(tmpfolder)
	
		shutil.copy(gmxfitness.files_route+gmxfitness.protein_file,tmpfolder)
		shutil.copy(gmxfitness.files_route+gmxfitness.ligand_name+'.itp',tmpfolder)
		shutil.copy(gmxfitness.files_route+gmxfitness.ligand_name+'.pdb',tmpfolder)
		shutil.copy(gmxfitness.files_route+gmxfitness.em_file,tmpfolder)
		shutil.copy(gmxfitness.files_route+gmxfitness.protein_ligand_box_file,tmpfolder)

		os.chdir(tmpfolder)  		

		try:
			n = len(gmxfitness.CMD)
			i = 0
			while(i<n):								
				p = subprocess.Popen(shlex.split(gmxfitness.CMD[i]), universal_newlines=True,
			                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				out,err = p.communicate()
				if i == 0:
					gmxfitness.process_topology()			
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
			if err is not None:		
				print err
			print "Error: %s" % e

		return energy

#Ejemplo
print gmxfitness.calculate_fitness('HiloPrueba')