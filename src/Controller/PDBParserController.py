from Tkinter import * 
from tkFileDialog import *
from Model.PDBParser import PDBParser

class PDBParserController:
	def __init__(self):
		self.pdbparser = PDBParser();			

	def openFile(self):		
		root = Tk()
		root.withdraw()
		return askopenfilename(**self.pdbparser.OPTIONS)		

	def readFile(self,protein_file_path):	
		file = open(protein_file_path,'r')
		while 1:
			line = file.readline()
			if not line: break
			items = line.split()

			if items[0] in self.pdbparser.REQUIREDRECORDTYPES:			
				self.pdbparser.record_type_found[self.pdbparser.REQUIREDRECORDTYPES.index(items[0])] = 1

			if items[0] == self.pdbparser.COORDINATERECORDTYPES[2]: #Helix
				self.getHelixSheetCoordenates(items,1)				
			elif items[0] == self.pdbparser.COORDINATERECORDTYPES[3]: #Sheet							
				self.getHelixSheetCoordenates(items,2)

			if items[0] == self.pdbparser.COORDINATERECORDTYPES[0]: #Atom 					
				if int(items[1]) in self.pdbparser.atomCoordenates:
					print "hola" # guardar las coordenadas del atomo si es que se encuentra en helix o sheet

		if 0 in self.pdbparser.record_type_found:
			IsValid = 0

	def validateFile(self):
		if self.pdbparser.IsValid:
			print "Valido"
		else:
			print "Invalido"			
	def getHelixSheetCoordenates(self,items,helixOrSheet):#1 para helix 2 para sheets
		try:
			if(helixOrSheet == 1):
				self.pdbparser.atomCoordenates.add(int(items[5]))
				self.pdbparser.atomCoordenates.add(int(items[8]))				
			elif(helixOrSheet == 2):				
				self.pdbparser.atomCoordenates.add(int(items[6]))
				self.pdbparser.atomCoordenates.add(int(items[9]))			
		except ValueError:
			self.IsValid = 0

