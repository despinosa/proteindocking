from Tkinter import * 
from tkFileDialog import *

class PDBParser:	
	OPTIONS = {}
	OPTIONS['defaultextension'] = '.pdb'
	OPTIONS['filetypes'] = [('PDB files','.pdb'),('ENT files','.ent'),('All files', '.*')]
	REQUIREDRECORDTYPES = ['HEADER',
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
				   'HELIX',
				   'SHEET']
	COORDINATERECORDTYPES = ['ATOM',
							'HETATM',
							'HELIX',
							'SHEET']
	record_type_found = [0]*len(REQUIREDRECORDTYPES)
	atomCoordenates = set()
	IsValid = 1    

	def openFile(self):		
		root = Tk()
		root.withdraw()
		return askopenfilename(**self.OPTIONS)		

	def readFile(self,protein_file_path):	
		file = open(protein_file_path,'r')
		while 1:
			line = file.readline()
			if not line: break
			items = line.split()

			if items[0] in self.REQUIREDRECORDTYPES:			
				self.record_type_found[self.REQUIREDRECORDTYPES.index(items[0])] = 1

			if items[0] == self.COORDINATERECORDTYPES[2]: #Helix
				self.getHelixSheetCoordenates(items,1)				
			elif items[0] == self.COORDINATERECORDTYPES[3]: #Sheet							
				self.getHelixSheetCoordenates(items,2)

			if items[0] == self.COORDINATERECORDTYPES[0]: #Atom 					
				if int(items[1]) in self.atomCoordenates:
					print "hola" # guardar las coordenadas del atomo si es que se encuentra en helix o sheet

		if 0 in self.record_type_found:
			IsValid = 0

	def validateFile(self):
		if self.IsValid:
			print "Valido"
		else:
			print "Invalido"			
	def getHelixSheetCoordenates(self,items,helixOrSheet):#1 para helix 2 para sheets
		try:
			if(helixOrSheet == 1):
				self.atomCoordenates.add(int(items[5]))
				self.atomCoordenates.add(int(items[8]))				
			elif(helixOrSheet == 2):				
				self.atomCoordenates.add(int(items[6]))
				self.atomCoordenates.add(int(items[9]))			
		except ValueError:
			self.IsValid = 0

pdbparser = PDBParser()
pdbfile = pdbparser.openFile()
pdbparser.readFile(pdbfile)
pdbparser.validateFile()
