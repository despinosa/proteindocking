import subprocess
import shutil
import os
import sys

class gmxfitness():

	def __init__(self,thread_name,protein_file,ligand_name,protein_ligand_box_file,em_file,forcefield):
		self.tmpfolder = thread_name
		self.protein_file = protein_file
		self.ligand_name = ligand_name
		self.protein_ligand_box_file = protein_ligand_box_file
		self.em = em_file		
		self.forcefield = forcefield		
		self.files_route = '../tmp/'			
		self.cmd = "gmx pdb2gmx -ignh -f "+ self.protein_file + " -ff " + forcefield + " -water none"		
		self.cmdgrompp = "gmx grompp -v -f em.mdp -c conf_with_ligand10x10x10.pdb -o em.tpr -p topol_with_ligand.top"
		self.cmdrun = "gmx mdrun -v -s em.tpr"

	
	def nspaces(self,n):
		cad = ''
		while(n):
			cad += ' '
			n -= 1
		return cad
	
	def ioFile(self):
		with open('topol_with_ligand.top','r') as in_file:
			buf = in_file.readlines()
		with open('topol_with_ligand.top','w') as out_file:	
			allowed = 0		
			newfile = ''
			for line in buf:			
				if '#include "charmm27.ff/forcefield.itp"' in line:
					line += '#include "'+self.ligand_name+'.itp"\n'				
				if '[ molecules ]' in line:
					allowed = 1
				if 'Protein_chain_A' in line and allowed :
					pos = line.find('1')
					if pos == -1:
						print 'Error'
						break
					line += (self.ligand_name+self.nspaces(pos-len(self.ligand_name)))+'1'		
				newfile += line
			out_file.write(newfile)
	
	def calculate_fitness(self):
		energy = ''		

		if os.path.exists(self.tmpfolder):
			shutil.rmtree(self.tmpfolder)

		os.mkdir(self.tmpfolder)
	
		shutil.copy(self.files_route+self.protein_file,self.tmpfolder)
		shutil.copy(self.files_route+self.ligand_name+'.itp',self.tmpfolder)
		shutil.copy(self.files_route+self.ligand_name+'.pdb',self.tmpfolder)
		shutil.copy(self.files_route+self.em,self.tmpfolder)
		shutil.copy(self.files_route+self.protein_ligand_box_file,self.tmpfolder)

		os.chdir(self.tmpfolder)  		

		try:
			p = subprocess.Popen(self.cmd, universal_newlines=True,
			                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out,err = p.communicate()

			if p.returncode:
				raise Exception(p.returncode)
			
			shutil.copyfile('topol.top','topol_with_ligand.top')
			self.ioFile()			

			p = subprocess.Popen(self.cmdgrompp, universal_newlines=True,
			                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out,err = p.communicate()

			if p.returncode:
				raise Exception(p.returncode)

			p = subprocess.Popen(self.cmdrun, universal_newlines=True,
			                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out,err = p.communicate()

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
			
			if os.path.exists(self.tmpfolder):
				shutil.rmtree(self.tmpfolder)

			if p.returncode:
				raise Exception(p.returncode)				
		except:
			e = sys.exc_info()[1]
			print err
			print "Error: %s" % e

		return energy

#Ejemplo
# aux = gmxfitness('HiloPrueba','1ubq.pdb','zinc_11909586','conf_with_ligand10x10x10.pdb','em.mdp','charmm27')
# print aux.calculate_fitness()