import subprocess
import shutil
import os

def nspaces(n):
	cad = ''
	while(n):
		cad += ' '
		n -= 1
	return cad

def ioFile(ligand_name):
	with open('topol_with_ligand.top','r') as in_file:
		buf = in_file.readlines()
	with open('topol_with_ligand.top','w') as out_file:	
		allowed = 0		
		newfile = ''
		for line in buf:			
			if '#include "charmm27.ff/forcefield.itp"' in line:
				line += '#include "'+ligand_name+'.itp"\n'				
			if '[ molecules ]' in line:
				allowed = 1
			if 'Protein_chain_A' in line and allowed :
				pos = line.find('1')
				if pos == -1:
					print 'Error'
					break
				line += (ligand_name+nspaces(pos-len(ligand_name)))+'1'		
			newfile += line
		out_file.write(newfile)

if os.path.exists('calc'):
	shutil.rmtree('calc')

os.mkdir('calc')

tmpfolder = 'calc'
file_name = "1ubq.pdb"
ligand_name = "zinc_11909586"
em = 'em.mdp'
shutil.copy('../tmp/'+file_name,tmpfolder)
shutil.copy('../tmp/'+ligand_name+'.itp',tmpfolder)
shutil.copy('../tmp/'+ligand_name+'.pdb',tmpfolder)
shutil.copy('../tmp/'+em,tmpfolder)

os.chdir('calc')  

forcefield = "charmm27"
cmd = "gmx pdb2gmx -ignh -f "+ file_name + " -ff " + forcefield + " -water none"
cmdbox = "gmx insert-molecules -f conf.gro -ci "+ligand_name+".pdb -o conf_with_ligand10x10x10.pdb -box 10 10 10 -nmol 1"
cmdgrompp = "gmx grompp -v -f em.mdp -c conf_with_ligand10x10x10.pdb -o em.tpr -p topol_with_ligand.top"
cmdrun = "gmx mdrun -v -s em.tpr"

try:

	p = subprocess.Popen(cmd, universal_newlines=True,
	                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out,err = p.communicate()

	if p.returncode:
		raise Exception(p.returncode)
	
	shutil.copyfile('topol.top','topol_with_ligand.top')
	ioFile(ligand_name)

	p = subprocess.Popen(cmdbox, universal_newlines=True,
	                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out,err = p.communicate()

	if p.returncode:
		raise Exception(p.returncode)

	p = subprocess.Popen(cmdgrompp, universal_newlines=True,
	                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out,err = p.communicate()

	if p.returncode:
		raise Exception(p.returncode)

	p = subprocess.Popen(cmdrun, universal_newlines=True,
	                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out,err = p.communicate()

	#obtiene la energia
	energy = ''
	pos = err.find('Epot=')
	if pos != -1:
		i = pos + 5
		while(1):
			if err[i] == 'F':
				break
			if err[i] != ' ':
				energy += err[i]		
			i += 1

	print energy
	os.chdir('..')
	
	if os.path.exists('calc'):
		shutil.rmtree('calc')

	if p.returncode:
		raise Exception(p.returncode)

except:
	e = sys.exc_info()[1]
	print "Error: %s" % e

