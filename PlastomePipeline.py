import sys
import os
import shutil
import gzip
import subprocess
import tkinter as tk
from tkinter import filedialog

'''
Housekeeping
'''

current_directory: str = os.getcwd()

'''
AT CLEANING ALGORITHM
'''

def ATCLEANER(InputFolder: str, CutoffMin: float, CutoffMax: float, MinRead: float):


    print("------- STARTING ATCLEANER ---------")

    Input_Folder: str = InputFolder # Name of the folder holding all the raw reads
    Cutoff_Percentage_From: float  = CutoffMin # Discard entries with an AT percentage below this ammount
    Cutoff_Percentage_To: float  = CutoffMax # Discard entries with an AT percentage above this ammount
    Minimuim_Read_Length: float = MinRead # Discard entries with a read length below this value
    Output_Folder: str = current_directory + "/ATCleanerOutput" # name of output directory

    if os.path.isdir(Output_Folder):
        for f in os.listdir(Output_Folder):
            os.remove(Output_Folder + "/" + f)
        os.removedirs(Output_Folder)

        os.makedirs(Output_Folder)
    else:
        os.makedirs(Output_Folder)

    print("\nAT cleaner V2.1 \nMade by Leonardo Cherin, UCL\n\n")
    print("Copyright 2025 Leonardo Cherin. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.\nThis program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.\nYou should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. \n\n\n")

    '''
    FASTQ format:

    @ = identifier for read
    ATAGCGA = sequence, all in capitals
    + = identifier for quality score initiation
    "£$"$%"$^ = representation of quality for each base called

    repeat until end
    '''

    #Input_File: str = input("Name of input file:\n")
    #Input_File: str = "example.fastq"
    #Cutoff_Percentage: float = float(input("What AT percentage cutoff should be used?\neg \"40\" indicates all sequences above 40% AT content will be retained, the rest will be discarded.\n"))
    #Cutoff_Percentage: float = 50.00

    Total_Entries: int = 0 # How many total reads
    Accepted_Entries: int = 0 # How any chloroplast reads have we kept?
    print("------- STARTED CLEANING ---------")

    for Input_File in os.listdir(str(Input_Folder)):
        subprocess.run(["gzip", "-dk", Input_Folder+"/"+Input_File]) # Unzip the raw read files
        with open(InputFolder+"/"+Input_File[:-3]) as file: # Open the file
            print("Processing File: " + str(Input_File))
            raw_data: list[str] = file.readlines() # Read the file as a list of lines
            location: int = 0 # Keep track of which line in the file I am
            for line in raw_data: # For each line
                if line[0] == "@": # Does it start with an identifier of a new read?
                    Total_Entries += 1 # Increase total reads by 1
                    sequence: str = raw_data[location+1] # look at the next line, this is the sequence!
                    sequence_length: int = 0 # keep track of total length of sequence
                    at_content: int = 0 # Keep track of the number of Adenines or Thymines encountered
                    for base in sequence:
                        sequence_length += 1
                        if base == "A" or base == "T" or base == "U":
                            at_content += 1
                    at_percentage: float = (at_content/sequence_length) * 100
                    with open(Output_Folder+"/ATCLEANER_LOG_"+Input_File.rpartition("/")[2][0:-6] + ".txt", "a") as logfile:
                        logfile.writelines("entry number " + str(Total_Entries) + " AT percentage: " + str(at_percentage) + " Read Length: " + str(sequence_length) + "\n")
                    if at_percentage >= Cutoff_Percentage_From and at_percentage < Cutoff_Percentage_To and sequence_length >= Minimuim_Read_Length:
                        with open(Output_Folder + "/" + Input_File.rpartition("/")[2][0:-6] + "_AT-cleaned.fastq", "a") as output:
                            output.writelines(raw_data[location:location+4]) # Write the raw entry to a new file, appending to the end. If the file does not exist, it creates it.
                        with open(Output_Folder+"/ATCLEANER_LOG_"+Input_File.rpartition("/")[2][0:-6] + ".txt", "a") as logfile:
                            logfile.writelines("ACCEPTED\n")
                        Accepted_Entries += 1
                    else:
                        with open(Output_Folder+"/ATCLEANER_LOG_"+Input_File.rpartition("/")[2][0:-6] + ".txt", "a") as logfile:
                            logfile.writelines("REJECTED\n")
                location += 1 # Move to the next line
            os.remove(Input_Folder + "/" + Input_File[:-3]) # Delete the un-compressed file
        with open(Output_Folder+"/ATCLEANER_LOG_"+Input_File.rpartition("/")[2][0:-6] + ".txt", "a") as logfile:
            logfile.writelines("Total Reads: " + str(Total_Entries) + "\n")
            logfile.writelines("Accepted Reads: " + str(Accepted_Entries) + "\n")
            logfile.writelines("Percentage of total kept: " + str((Accepted_Entries/Total_Entries)*100) + "\n")
            logfile.writelines("Output file name: " + Output_Folder + Input_File[0:-6] + "_AT-cleaned.fastq" + "\n\n")


    print("------- FINISHED CLEANING ---------")
    print("------- COMPRESSING OUTPUT---------")
    subprocess.run(["gzip -r " + Output_Folder + "/*.fastq"], shell=True) # Re-compress the output files
    ATcleanerhasrun.set(True)

    print("------- ATCLEANER FINISHED ---------")


def PlastomeAlignment(ReferenceGenomePath: str):
    Output_Folder: str = current_directory + "/PlastomeAlignmentOutput" # name of output directory
    Input_Folder: str = current_directory + "/ATCleanerOutput" # name of output directory

    if os.path.isdir(Output_Folder):
        for f in os.listdir(Output_Folder):
            os.remove(Output_Folder + "/" + f)
        os.removedirs(Output_Folder)

        os.makedirs(Output_Folder)
    else:
        os.makedirs(Output_Folder)

    subprocess.run(["minimap2 -ax map-ont " + ReferenceGenomePath + " " + Input_Folder+ "/*.fastq.gz" + " > " + Output_Folder + "/alignment.sam"], shell=True)

    subprocess.run(["samtools", "view", "-@", "n", "-Sb", "-o", Output_Folder+"/alignment.bam", Output_Folder+"/alignment.sam"])

    subprocess.run(["samtools", "sort", "-O", "bam", "-o", Output_Folder+"/sorted_alignment.bam", Output_Folder+"/alignment.bam"])

    subprocess.run(["samtools", "index", Output_Folder+"/sorted_alignment.bam"])

    os.remove(Output_Folder+"/alignment.bam")
    os.remove(Output_Folder+"/alignment.sam")

def PlastomeAssemble(ReferenceGenomePath: str):
    Input_Folder: str = current_directory + "/ATCleanerOutput"
    Output_Folder: str = current_directory + "/PlastomeDeNovoAssemblyOutput"
    
    if os.path.isdir(Output_Folder):
        shutil.rmtree(Output_Folder)
        os.makedirs(Output_Folder)
        os.makedirs(Output_Folder+"/AlignmentCheck")
        os.makedirs(Output_Folder+"/PorechopCleaned")
    else:
        os.makedirs(Output_Folder)
        os.makedirs(Output_Folder+"/AlignmentCheck")
        os.makedirs(Output_Folder+"/PorechopCleaned")

    for file in os.listdir(Input_Folder):
        if file[-9:] == ".fastq.gz":
            subprocess.run(["porechop", "-i", Input_Folder+"/"+file, "-o", Output_Folder+"/PorechopCleaned/PORECHOP_"+file])
    #subprocess.run(["cat "+Output_Folder+"/PorechopCleaned/*.fastq.gz > "+ Output_Folder+"/concatenated_porechop_cleaned.fastq.gz"], shell=True)
    #subprocess.run(["flye --nano-hq " + Output_Folder + "/concatenated_porechop_cleaned.fastq.gz --out-dir " + Output_Folder + "/FlyeOutput --genome-size 0.0002g -i 3"], shell=True)
    shutil.rmtree(Output_Folder+"/PorechopCleaned")
    os.remove(Output_Folder+"/concatenated_porechop_cleaned.fastq.gz")
    
    if ReferenceGenomePath != "":
        subprocess.run(["minimap2 -ax map-ont " + ReferenceGenomePath + " " + Output_Folder+ "/FlyeOutput/assembly.fasta" + " > " + Output_Folder + "/AlignmentCheck/alignment.sam"], shell=True)
        subprocess.run(["samtools", "view", "-@", "n", "-Sb", "-o", Output_Folder+"/AlignmentCheck/alignment.bam", Output_Folder+"/AlignmentCheck/alignment.sam"])

        subprocess.run(["samtools", "sort", "-O", "bam", "-o", Output_Folder+"/AlignmentCheck/sorted_alignment.bam", Output_Folder+"/AlignmentCheck/alignment.bam"])

        subprocess.run(["samtools", "index", Output_Folder+"/AlignmentCheck/sorted_alignment.bam"])

        os.remove(Output_Folder+"/AlignmentCheck/alignment.sam")
        os.remove(Output_Folder+"/AlignmentCheck/alignment.bam")

'''
GUI STUFF
'''

# Requires chmod u+rx

def Clean():
    ATCLEANER(input_settings["inputdir"], float(ATmin.get()), float(ATmax.get()), float(minReadLength.get()))

def Alignment():
    print("------- STARTING ALIGNMENT ---------")
    if ATcleanerhasrun.get() == False:
        Clean()
    PlastomeAlignment(input_settings["genomepath"]) 
    print("------- FINISHED ALIGNMENT ---------") 

def Assembly():
    if ATcleanerhasrun.get() == False:
        Clean()
    print("-------- STARTING ASSEMBLY ---------") 
    PlastomeAssemble(input_settings["genomepath"])
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
ATcleanerhasrun = tk.BooleanVar()


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
minreadlengthlabel = tk.Label(root, text="minimum read length (bases):").pack()
MinReadLength = tk.Entry(root, textvariable=minReadLength)
MinReadLength.insert(0,"10000")
MinReadLength.pack()

referencegenomelabel = tk.Label(root, text="Select reference genome (Align only):").pack()
SelectInput = tk.Button(root, text="Reference Genome", command=GenomeLocation).pack()

startprocesseslabel = tk.Label(root, text="Run pipelines:").pack()
LaunchClean = tk.Button(root, text="ATCleaner", command=Clean).pack()
CleanCheck = tk.Checkbutton(root, text="Has AT cleaner run?", variable=ATcleanerhasrun).pack()
LaunchAlign = tk.Button(root, text="Alignment Pipeline", command=Alignment).pack()
LaunchAssemble = tk.Button(root, text="Assembly Pipeline", command=Assembly).pack()

root.mainloop()
