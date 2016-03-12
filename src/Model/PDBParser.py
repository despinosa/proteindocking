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
	def __init__(self):		
		self.atomCoordenates = set()
		self.record_type_found = [0]*len(self.REQUIREDRECORDTYPES)	
		self.IsValid = 1
