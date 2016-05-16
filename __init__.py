import os
import sys
from src.gui.GUI import GUI

def __init__(self): 	
	path = os.path.dirname(__file__)
	sys.path.append(path)
	self.menuBar.addmenuitem('Plugin', 'command',
	                         'Protein docking',
	                         label = 'Protein docking',
	                         command = lambda s=self : GUI(s))