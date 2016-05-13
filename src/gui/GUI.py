from Tkinter import * 
from tkFileDialog import askopenfilename

class GUI(Toplevel):
    def __init__(self,master):                    
        Toplevel.__init__(self,master.root,height=250,width=650)
        self.resizable(width=FALSE,height=FALSE)        
        self.wm_title('Protein docking')
        self.center()
        #Labels
        self.l_protein_pdb = Label(self, text="Protein PDB file: ",padx=10,pady=10)
        self.l_ligand_pdb = Label(self, text="Ligand PDB file: ",padx=10,pady=10)
        self.l_ligand_itp = Label(self, text="Ligand ITP file: ",padx=10,pady=10)
        self.l_cavities = Label(self, text="Cavities PDB file: ",padx=10,pady=10)
        self.l_protein_pdb.pack(side="left")
        self.l_ligand_pdb.pack(side="left")
        self.l_ligand_itp.pack(side="left")
        self.l_cavities.pack(side="left")
        self.l_protein_pdb.place(x=10,y=10)        
        self.l_ligand_pdb.place(x=10,y=50)        
        self.l_ligand_itp.place(x=10,y=100)        
        self.l_cavities.place(x=10,y=150)  
        #Buttons
        self.b_protein_pdb = Button(self,text="Choose file...",name="b_protein_pdb")
        self.b_protein_pdb.bind("<Button-1>",self.openFile)
        self.b_ligand_pdb = Button(self,text="Choose file...",name="b_ligand_pdb")
        self.b_ligand_pdb.bind("<Button-1>",self.openFile)
        self.b_ligand_itp = Button(self,text="Choose file...",name="b_ligand_itp")
        self.b_ligand_itp.bind("<Button-1>",self.openFile)
        self.b_cavities = Button(self,text="Choose file...",name="b_cavities")
        self.b_cavities.bind("<Button-1>",self.openFile)
        self.b_protein_pdb.pack(side="right")
        self.b_ligand_pdb.pack(side="right")
        self.b_ligand_itp.pack(side="right")
        self.b_cavities.pack(side="right")
        self.b_protein_pdb.place(x=550,y=12)
        self.b_ligand_pdb.place(x=550,y=52)
        self.b_ligand_itp.place(x=550,y=102)
        self.b_cavities.place(x=550,y=152)
        #Textboxes
        self.txt_protein_pdb = Text(self,height=1, width=70,stat=DISABLED)
        self.txt_ligand_pdb = Text(self,height=1, width=70,state=DISABLED)
        self.txt_ligand_itp = Text(self,height=1, width=70,state=DISABLED)
        self.txt_cavities = Text(self,height=1, width=70,state=DISABLED)
        self.txt_protein_pdb.pack(side="top")
        self.txt_ligand_pdb.pack(side="top")
        self.txt_ligand_itp.pack(side="top")
        self.txt_cavities.pack(side="top")
        self.txt_protein_pdb.place(x=110,y=18)
        self.txt_ligand_pdb.place(x=110,y=58)
        self.txt_ligand_itp.place(x=110,y=108)
        self.txt_cavities.place(x=110,y=158)
        #Docking button
        self.run_docking = Button(self,text="Start docking",name="run_docking",state=DISABLED)
        self.run_docking.pack(side="bottom")      
        self.run_docking.place(x=self.size[0]/2 - 50,y=self.size[1] - 50)         

        self.lift()

    def openFile(self,event):
        OPTIONS = {}
        OPTIONS['defaultextension'] = '.pdb'                 
        OPTIONS['filetypes'] = ([('ITP files','.itp')] if self.b_ligand_itp == event.widget else [('PDB files','.pdb')])
        file_name = askopenfilename(**OPTIONS)                
        if(self.b_protein_pdb == event.widget):
            self.txt_protein_pdb["state"] = "normal"
            self.txt_protein_pdb.delete(1.0,END)
            self.txt_protein_pdb.insert(INSERT,file_name)
            self.txt_protein_pdb["state"] = "disabled"  
        elif(self.b_ligand_pdb == event.widget):
            self.txt_ligand_pdb["state"] = "normal"
            self.txt_ligand_pdb.delete(1.0,END)
            self.txt_ligand_pdb.insert(INSERT,file_name)
            self.txt_ligand_pdb["state"] = "disabled"  
        elif(self.b_ligand_itp == event.widget):                        
            self.txt_ligand_itp["state"] = "normal"
            self.txt_ligand_itp.delete(1.0,END)
            self.txt_ligand_itp.insert(INSERT,file_name)
            self.txt_ligand_itp["state"] = "disabled"  
        elif(self.b_cavities == event.widget):
            self.txt_cavities["state"] = "normal"
            self.txt_cavities.delete(1.0,END)
            self.txt_cavities.insert(INSERT,file_name)                
            self.txt_cavities["state"] = "disabled"  

    #center frame
    def center(self):
        self.update_idletasks()
        w = self.winfo_screenwidth()
        h = self.winfo_screenheight()
        self.size = tuple(int(_) for _ in self.geometry().split('+')[0].split('x'))
        self.x = w/2 - self.size[0]/2
        self.y = h/2 - self.size[1]/2
        self.geometry("%dx%d+%d+%d" % (self.size + (self.x, self.y)))

