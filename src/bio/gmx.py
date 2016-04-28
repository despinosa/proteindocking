import subprocess
import shutil
import os
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
			
	dockedpair = 'dockedpair.pdb'
	em_file = 'em.mdp'			
	forcefield = 'charmm27'
	topol_with_ligand_file = 'topol_with_ligand.top' 		

	@staticmethod
	def nspaces(n):
		cad = ''
		while(n):
			cad += ' '
			n -= 1
		return cad
	@staticmethod
	def add_molecule(line,ligand_name):
		pos = line.find('1')
		if pos == -1:
			raise Exception('Error')
		return (ligand_name+gmx.nspaces(pos-len(ligand_name)))+'1'
	@staticmethod	
	def ioFile(dp_object):		
		with open(gmx.topol_with_ligand_file,'r') as in_file:
			buf = in_file.readlines()
		with open(gmx.topol_with_ligand_file,'w') as out_file:	
			allowed = 0		
			newfile = ''
			for i,line in enumerate(buf):			
				if '#include "'+gmx.forcefield+'.ff/forcefield.itp"' in line:
					line += '#include "'+dp_object.ligand_name+'.itp"\n'				
				if '[ molecules ]' in line:
					allowed = 1
				if 'Protein' in line and allowed:					
					if i == len(buf) - 1:											
						line += gmx.add_molecule(line,dp_object.ligand_name)
					elif '[' or ']' in line:						
						line += gmx.add_molecule(line,dp_object.ligand_name)
				newfile += line
			out_file.write(newfile)		

	@staticmethod
	def process_topology():
		os.chdir(os.path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES))
		shutil.copyfile('topol.top',gmx.topol_with_ligand_file)
		gmx.ioFile()

	@staticmethod
	def generate_protein_topology(dp_object):						
		os.chdir(os.path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES))
		cmd_protein_topology = "gmx pdb2gmx -ignh -f {} -ff {} -water none".format(dp_object.protein_file,gmx.forcefield)
		try:
			p = subprocess.Popen(shlex.split(cmd_protein_topology), universal_newlines=True,
				                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out,err = p.communicate()
		except Exception:
			e = sys.exc_info()[1]
			print "Error: %s" % e
	# @staticmethod
	# def protein_ligand_box():
	# 	thread_name = 'Hilo0'
	# 	#thread_name = current_thread().name		
	# 	os.chdir(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.TMP,thread_name,gmx.GMX_FILES))		
	# 	shutil.copyfile(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.TMP,thread_name,gmx.dockedpair),gmx.dockedpair)		
	# 	cmd_protein_ligand_box = "gmx insert-molecules -f conf.gro -ci {} -o {} -box {} {} {} -nmol 1".format(gmx.dockedpair,gmx.protein_ligand_box_file,gmx.w,gmx.h,gmx.l)		
	# 	try:
	# 		p = subprocess.Popen(shlex.split(cmd_protein_ligand_box), universal_newlines=True,
	# 			                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# 		out,err = p.communicate()
	# 		print out
	# 		print err
	# 	except Exception:
	# 		e = sys.exc_info()[1]				
	# 		print "Error: %s" % e

	@staticmethod
	def add_hydrogens(dp_object):		
		os.chdir(os.path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES))
		cmd_hydrogens = ["gmx grompp -f {} -o em_aux.tpr -c conf.gro".format(gmx.em_file),
					     "gmx trjconv -f conf.gro -o {} -s em_aux.tpr".format(dp_object.protein_file)]
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
				p.stdin.close()						
				i += 1
		except Exception:
			e = sys.exc_info()[1]				
			print "Error: %s" % e
	@staticmethod
	def make_temp_folders():
            if os.path.exists(os.path.join(gmx.TEMPDIR,gmx.ROOT)):
                shutil.rmtree(os.path.join(gmx.TEMPDIR,gmx.ROOT))
            os.mkdir(os.path.join(gmx.TEMPDIR,gmx.ROOT))
            os.chdir(os.path.join(gmx.TEMPDIR,gmx.ROOT))
            os.mkdir(gmx.FILES)
            os.mkdir(gmx.TMP)

	@staticmethod
	def process_folders(dp_object):			
		thread_name = current_thread().name
		#thread_name = 'Hilo0'

		os.chdir(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.TMP))		

		if not os.path.exists(thread_name):			
			os.mkdir(thread_name)

		os.chdir(thread_name)

		if os.path.exists(gmx.GMX_FILES):
			shutil.rmtree(gmx.GMX_FILES)

		os.mkdir(gmx.GMX_FILES)
	
		shutil.copy(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.FILES, dp_object.protein_file), gmx.GMX_FILES)		
		shutil.copy(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.FILES, '{}.itp'.format(dp_object.ligand_name)), gmx.GMX_FILES)
		shutil.copy(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.FILES, '{}.pdb'.format(dp_object.ligand_name)), gmx.GMX_FILES)
		shutil.copy(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.FILES, gmx.em_file), gmx.GMX_FILES)
		shutil.copy(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.FILES, gmx.topol_with_ligand_file), gmx.GMX_FILES)				
		shutil.copy(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.FILES, 'conf.gro'), gmx.GMX_FILES)				


		os.chdir(gmx.GMX_FILES)

	@staticmethod
	def calculate_fitness():
		energy = ''			
		#thread_name = current_thread().name
		thread_name = 'Hilo0'
		os.chdir(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.TMP,thread_name,gmx.GMX_FILES))	
		shutil.copyfile(os.path.join(gmx.TEMPDIR, gmx.ROOT, gmx.TMP,thread_name,gmx.dockedpair),gmx.dockedpair)	
		cmd_fitness = ["gmx grompp -v -f {} -c {} -o em.tpr -p {}".format(gmx.em_file,gmx.dockedpair,gmx.topol_with_ligand_file),
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
		return energy

if __name__ == '__main__':
	#script_route = os.getcwd()
	#gmx.make_temp_folders()
	# gmx.generate_protein_topology()	
	# os.chdir(os.path.join(gmx.TEMPDIR,gmx.ROOT,gmx.FILES))
	# os.rename(gmx.protein_file,'{}_'.format(gmx.protein_file))
	# gmx.add_hydrogens()
	# gmx.generate_protein_topology()
	#gmx.process_topology()
	#gmx.process_folders()
	#gmx.protein_ligand_box()
	print float(gmx.calculate_fitness())


