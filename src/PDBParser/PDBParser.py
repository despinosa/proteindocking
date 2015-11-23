from Tkinter import * 
from tkFileDialog import *

class PDBParser:	
	options = {}
	options['defaultextension'] = '.pdb'
	options['filetypes'] = [('PDB files','.pdb'),('All files', '.*')]
	record_type_found = [0]*21
	RECORDTYPES = ['HEADER',
				   'TITLE',
				   'COMPND',
				   'SOURCE',
				   'KEYWDS',
				   'EXPDTA',
				   'AUTHOR',
				   'REVDAT',
				   #'REMARK 2',
				   #'REMARK 3',
				   'REMARK',
				   'SEQRES',
				   'CRYST1',
				   'ORIGX1',																																																																																																								
				   'ORIGX2',
				   'ORIGX3',
				   'SCALE1',
				   'SCALE2',
				   'SCALE3',
				   'MASTER',
				   'END',
				   'HELIX', # ???
				   'SHEET'] # ???
    
	@staticmethod
	def openFile():		
		root = Tk()
		root.withdraw()
		return askopenfilename(**PDBParser.options)		
	@staticmethod
	def readFile(protein_file_path):
		file = open(protein_file_path,'r')
		while 1:
			line = file.readline()
			if not line: break
			items = line.split()

			if items[0] in PDBParser.RECORDTYPES:
				print items[0]				
				PDBParser.record_type_found[PDBParser.RECORDTYPES.index(items[0])] = 1			
	@staticmethod 
	def validateFile():
		print PDBParser.record_type_found
		if 0 in PDBParser.record_type_found:
			print "Invalido"
		else:
			print "Valido"

pdbfile = PDBParser.openFile()
PDBParser.readFile(pdbfile)
PDBParser.validateFile()
