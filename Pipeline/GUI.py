import sys
import subprocess
import tkinter as tk
from tkinter import filedialog

# Requires chmod u+rx

def Alignment():
    print("------- STARTING ALIGNMENT ---------")
    subprocess.check_call(["./PlastomeAlignPipeline.sh", input_settings["inputdir"], ATmin.get(), ATmax.get(), minReadLength.get(), input_settings["genomepath"]])
    print("------- FINISHED ALIGNMENT ---------") 

def Assembly():
    print("-------- STARTING ASSEMBLY ---------") 
    subprocess.check_call(["./PlastomeAssemblePipeline.sh", input_settings["inputdir"], ATmin.get(), ATmax.get(), minReadLength.get()])
    print("-------- FINISHED ASSEMBLY ---------") 
def OpenInputFolder():
    InputFolder = filedialog.askdirectory()
    input_settings["inputdir"] = InputFolder
def GenomeLocation():
    ReferenceGenomeLocation = filedialog.askopenfilename(filetypes=[("FASTA","*.fasta *.fa")])
    input_settings["genomepath"] = ReferenceGenomeLocation
    
root = tk.Tk()

input_settings: dict = {
    "inputdir" : "",
    "genomepath" : "",
}

ATmin = tk.StringVar()
ATmax = tk.StringVar()
minReadLength = tk.StringVar()


root.title("Plastome Pipeline")

mainframe = tk.Frame(root)

TitleLabel = tk.Label(root, text="Plastome Pipeline").pack()

inputfolderlabel = tk.Label(root, text="Select Input Folder:").pack()
SelectInput = tk.Button(root, text="Input Folder", command=OpenInputFolder).pack()

minATlabel = tk.Label(root, text="AT % minimum:").pack()
MinAT = tk.Entry(root, textvariable=ATmin)
MinAT.insert(0,"57")
MinAT.pack()
maxATlabel = tk.Label(root, text="AT % maximum:").pack()
MaxAT = tk.Entry(root, textvariable=ATmax)
MaxAT.insert(0,"100")
MaxAT.pack()
minreadlengthlabel = tk.Label(root, text="minimum read length:").pack()
MinReadLength = tk.Entry(root, textvariable=minReadLength)
MinReadLength.insert(0,"10000")
MinReadLength.pack()

referencegenomelabel = tk.Label(root, text="Select reference genome (Align only):").pack()
SelectInput = tk.Button(root, text="Reference Genome", command=GenomeLocation).pack()

startprocesseslabel = tk.Label(root, text="Run pipelines:").pack()
LaunchAlign = tk.Button(root, text="Alignment Pipeline", command=Alignment).pack()
LaunchAssemble = tk.Button(root, text="Assembly Pipeline", command=Assembly).pack()

root.mainloop()
