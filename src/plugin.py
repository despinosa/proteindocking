from gui.GUI import GUI

def __init__(self): 	
	self.menuBar.addmenuitem('Plugin', 'command',
	                         'Protein docking',
	                         label = 'Protein docking',
	                         command = lambda s=self : GUI(s))