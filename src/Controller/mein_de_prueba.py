from PDBParserController import PDBParserController

aux = PDBParserController()
pdbfile = aux.openFile()
aux.readFile(pdbfile)
aux.validateFile()