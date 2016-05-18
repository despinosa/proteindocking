import os
import Queue
import time
import ttk, threading
import tkMessageBox
from datetime import datetime
from sys import stdout
from time import sleep    
from os import path 
from alpsmain import ALPSMain
from Tkinter import * 
from tkFileDialog import askopenfilename, askdirectory


class GUI(Toplevel):
    def __init__(self,master,files_path):                    
        Toplevel.__init__(self,master.root,height=400,width=700)        
        self.resizable(width=FALSE,height=FALSE)        
        self.wm_title('Protein docking')
        self.center()
        self.DIR_OPTIONS = {}        
        self.DIR_OPTIONS['mustexist'] = True        
        self.OPTIONS = {}
        self.OPTIONS['defaultextension'] = '.pdb'    
        self.FILES_PATH = files_path
        self.valid_protein = 0
        self.valid_ligand_pdb = 0
        self.valid_ligand_itp = 0
        self.valid_cavity = 0
        self.valid_output = 0                            
        self.valid_ff = 0
        self.forcefields = ['charmm27','gromos54a7_atb']        
        #Labels
        self.l_protein_pdb = Label(self, text="Protein PDB file: ",padx=10,pady=10)
        self.l_ligand_pdb = Label(self, text="Ligand PDB file: ",padx=10,pady=10)
        self.l_ligand_itp = Label(self, text="Ligand ITP file: ",padx=10,pady=10)
        self.l_cavities = Label(self, text="Cavities PDB file: ",padx=10,pady=10)
        self.l_output = Label(self, text="Output path: ",padx=10,pady=10)
        self.l_forcefield = Label(self, text="Forcefield: ",padx=10,pady=10)
        self.l_protein_pdb.pack(side="left")
        self.l_ligand_pdb.pack(side="left")
        self.l_ligand_itp.pack(side="left")
        self.l_cavities.pack(side="left")
        self.l_output.pack(side="left")
        self.l_forcefield.pack(side="left")
        self.l_protein_pdb.place(x=10,y=30)        
        self.l_ligand_pdb.place(x=10,y=70)        
        self.l_ligand_itp.place(x=10,y=120)        
        self.l_cavities.place(x=10,y=170)  
        self.l_output.place(x=10,y=220)
        self.l_forcefield.place(x=10,y=270)
        #Buttons
        self.b_protein_pdb = Button(self,text="Choose file...",name="b_protein_pdb")
        self.b_protein_pdb.bind("<Button-1>",self.openFile)
        self.b_ligand_pdb = Button(self,text="Choose file...",name="b_ligand_pdb")
        self.b_ligand_pdb.bind("<Button-1>",self.openFile)
        self.b_ligand_itp = Button(self,text="Choose file...",name="b_ligand_itp")
        self.b_ligand_itp.bind("<Button-1>",self.openFile)
        self.b_cavities = Button(self,text="Choose file...",name="b_cavities")
        self.b_cavities.bind("<Button-1>",self.openFile)
        self.b_output = Button(self,text="Choose folder...",name="b_output")
        self.b_output.bind("<Button-1>",self.getDir)
        self.b_protein_pdb.pack(side="right")
        self.b_ligand_pdb.pack(side="right")
        self.b_ligand_itp.pack(side="right")
        self.b_cavities.pack(side="right")
        self.b_output.pack(side="right")
        self.b_protein_pdb.place(x=550,y=32)
        self.b_ligand_pdb.place(x=550,y=72)
        self.b_ligand_itp.place(x=550,y=122)
        self.b_cavities.place(x=550,y=172)
        self.b_output.place(x=550,y=222)
        #Textboxes
        self.txt_protein_pdb = Text(self,height=1, width=70,stat=DISABLED)
        self.txt_ligand_pdb = Text(self,height=1, width=70,state=DISABLED)
        self.txt_ligand_itp = Text(self,height=1, width=70,state=DISABLED)
        self.txt_cavities = Text(self,height=1, width=70,state=DISABLED)
        self.txt_output = Text(self,height=1,width=70,state=DISABLED)
        self.txt_protein_pdb.pack(side="top")
        self.txt_ligand_pdb.pack(side="top")
        self.txt_ligand_itp.pack(side="top")
        self.txt_cavities.pack(side="top")
        self.txt_output.pack(side="top")
        self.txt_protein_pdb.place(x=110,y=38)
        self.txt_ligand_pdb.place(x=110,y=78)
        self.txt_ligand_itp.place(x=110,y=128)
        self.txt_cavities.place(x=110,y=178)
        self.txt_output.place(x=110,y=228)
        #Combobox
        self.cb_forcefields = ttk.Combobox(self,state="readonly")
        self.cb_forcefields.bind("<<ComboboxSelected>>",self.validate_ff)
        self.cb_forcefields["values"] = self.forcefields        
        self.cb_forcefields.pack(side="top")
        self.cb_forcefields.place(x=110,y=278)
        #Docking button
        self.run_docking = Button(self,text="Start docking",name="run_docking",state=DISABLED)        
        self.run_docking.pack(side="bottom")      
        self.run_docking.place(x=self.size[0]/2 - 50,y=self.size[1] - 50)                 

        self.lift()
    #Validations        
    def validate(self):                
        if(self.valid_protein and self.valid_ligand_pdb and self.valid_ligand_itp and self.valid_cavity and self.valid_output and self.valid_ff):            
            self.run_docking.bind("<Button-1>",self.start_docking)
            self.run_docking.config(state="normal")        

    def validate_ff(self,event):
        for i,ff in enumerate(self.forcefields):
            if(self.cb_forcefields.get() == ff):
                self.forcefield = i        
        self.valid_ff = 1
        self.validate()
    #Events
    def openFile(self,event):                
        self.OPTIONS['filetypes'] = ([('ITP files','.itp')] if self.b_ligand_itp == event.widget else [('PDB files','.pdb')])
        if(len(self.txt_protein_pdb.get(1.0,END)) > 0):            
            self.OPTIONS['initialdir'] = os.path.dirname(self.txt_protein_pdb.get(1.0,END))
        file_name = askopenfilename(**self.OPTIONS)                
        if(self.b_protein_pdb == event.widget):
            self.txt_protein_pdb["state"] = "normal"
            self.txt_protein_pdb.delete(1.0,END)
            self.txt_protein_pdb.insert(INSERT,file_name)
            self.txt_protein_pdb["state"] = "disabled"  
            self.protein_path = file_name
            self.valid_protein = 1
        elif(self.b_ligand_pdb == event.widget):
            self.txt_ligand_pdb["state"] = "normal"
            self.txt_ligand_pdb.delete(1.0,END)
            self.txt_ligand_pdb.insert(INSERT,file_name)
            self.txt_ligand_pdb["state"] = "disabled"  
            self.ligand_path = file_name
            self.valid_ligand_pdb = 1
        elif(self.b_ligand_itp == event.widget):                        
            self.txt_ligand_itp["state"] = "normal"
            self.txt_ligand_itp.delete(1.0,END)
            self.txt_ligand_itp.insert(INSERT,file_name)
            self.txt_ligand_itp["state"] = "disabled"  
            self.itp_path = file_name
            self.valid_ligand_itp = 1
        elif(self.b_cavities == event.widget):
            self.txt_cavities["state"] = "normal"
            self.txt_cavities.delete(1.0,END)
            self.txt_cavities.insert(INSERT,file_name)                
            self.txt_cavities["state"] = "disabled"
            self.cavities_path = file_name
            self.valid_cavity = 1
        self.validate()

    def getDir(self,event):        
        output_path = askdirectory(**self.DIR_OPTIONS)                
        if(self.b_output == event.widget):
            self.txt_output["state"] = "normal"
            self.txt_output.delete(1.0,END)
            self.txt_output.insert(INSERT,output_path)
            self.txt_output["state"] = "disabled"
            self.output = output_path
            self.valid_output = 1
        self.validate()

    def start_docking(self,event):                
        self.new_progress_bar()        
        self.queue = Queue.Queue()                  
        self.main = ALPSMain(self.queue,self.protein_path,self.ligand_path,self.itp_path,self.cavities_path,self.output,self.forcefield,self.FILES_PATH)
        self.main.start()
        self.after(100, self.process_queue)

    #Auxiliary functions
    def new_progress_bar(self):
        self.prog_bar = ttk.Progressbar(self, orient="horizontal",length=300, mode="determinate",maximum=100)
        self.prog_bar.pack(side=TOP)
        self.prog_bar["value"] = 0        
                   
    #center frame
    def center(self):
        self.update_idletasks()
        w = self.winfo_screenwidth()
        h = self.winfo_screenheight()
        self.size = tuple(int(_) for _ in self.geometry().split('+')[0].split('x'))
        self.x = w/2 - self.size[0]/2
        self.y = h/2 - self.size[1]/2
        self.geometry("%dx%d+%d+%d" % (self.size + (self.x, self.y)))        

    def process_queue(self):
        def check_errors():
            if not self.main.ex_queue.empty():
                e = self.main.ex_queue.get(0)
                tkMessageBox.showinfo("Protein docking", e)
                return True
            return False

        try:
            if check_errors(): return
            msg = self.queue.get(0)
            step_ = float(msg)                        
            self.prog_bar["value"] = step_
            self.after(100,self.process_queue)
        except ValueError:
            self.prog_bar.stop()
            tkMessageBox.showinfo("Protein docking",msg)
            self.destroy()                
        except Queue.Empty:
            self.after(100, self.process_queue)
