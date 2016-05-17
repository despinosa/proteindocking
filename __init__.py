import os
import sys
path_ = os.path.dirname(__file__)
sys.path.append(os.path.join(path_,'src'))	
from gui.GUI import GUI

def __init__(self): 		
	self.menuBar.addmenuitem('Plugin', 'command',
	                         'Protein docking',
	                         label = 'Protein docking',
	                         command = lambda s=self : GUI(s,path_))