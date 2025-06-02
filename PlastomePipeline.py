import sys
import os
import shutil
import gzip
import subprocess
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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

    ATandlengthDistribution = []

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
    Total_Accepted_Bases: int = 0 # how many bases have we got in our data?
    Accepted_Entries: int = 0 # How any chloroplast reads have we kept?
    print("------- STARTED CLEANING ---------")

    for Input_File in os.listdir(str(Input_Folder)):
        subprocess.run(["gzip", "-dk", Input_Folder+"/"+Input_File]) # Unzip the raw read files
        with open(Input_Folder+"/"+Input_File[:-3]) as file: # Open the file
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

                    if at_percentage >= Cutoff_Percentage_From and at_percentage < Cutoff_Percentage_To and sequence_length >= Minimuim_Read_Length: # Does the sequence satisfy the plastome constraints?
                        Total_Accepted_Bases += sequence_length
                        with open(Output_Folder + "/" + Input_File.rpartition("/")[2][0:-6] + "_AT-cleaned.fastq", "a") as output:
                            output.writelines(raw_data[location:location+4]) # Write the raw entry to a new file, appending to the end. If the file does not exist, it creates it.

                        with open(Output_Folder+"/ATCLEANER_LOG_"+Input_File.rpartition("/")[2][0:-6] + ".txt", "a") as logfile: # Create a log file entry for debugging purposes
                            logfile.writelines("ACCEPTED\n")

                        Accepted_Entries += 1 # We have accepted this sequence!

                    else: # Otherwise, the constraints were not satisfied
                        with open(Output_Folder+"/ATCLEANER_LOG_"+Input_File.rpartition("/")[2][0:-6] + ".txt", "a") as logfile: # Log the sequence as rejected
                            logfile.writelines("REJECTED\n")

                ATandlengthDistribution.append((at_percentage, sequence_length)) # Keep track of the length and AT content of all sequences for data visualisation

                location += 1 # Move to the next line
            os.remove(Input_Folder + "/" + Input_File[:-3]) # Delete the un-compressed file
        with open(Output_Folder+"/ATCLEANER_LOG_"+Input_File.rpartition("/")[2][0:-6] + ".txt", "a") as logfile:
            logfile.writelines("Total Reads: " + str(Total_Entries) + "\n")
            logfile.writelines("Accepted Reads: " + str(Accepted_Entries) + "\n")
            logfile.writelines("Percentage of total kept: " + str((Accepted_Entries/Total_Entries)*100) + "\n")
            logfile.writelines("Output file name: " + Output_Folder + Input_File[0:-6] + "_AT-cleaned.fastq" + "\n\n")

    print("------- FINISHED CLEANING ---------")
    print("------- COMPRESSING OUTPUT---------")
    subprocess.run(["gzip -r " + Output_Folder + "/*.fastq"], shell=True) # Compress the output files
    ATcleanerhasrun.set(True)

    if input_settings["genomepath"] != "":
        genomesize: int = 0
        with open(input_settings["genomepath"]) as file:
            for line in file.readlines():
                if line[0] != ">":
                    for character in line:
                        genomesize += 1

        coverage: float = Total_Accepted_Bases / genomesize
        Coverage.config(state=tk.NORMAL)
        Coverage.delete('1.0', tk.END)
        Coverage.insert(tk.END, str(coverage))
        Coverage.config(state=tk.DISABLED)
                        

    if os.path.isfile("SequenceDataDistribution.csv"):
        os.remove("SequenceDataDistribution.csv")
    outdata = pd.DataFrame(ATandlengthDistribution)
    outdata.columns = ["AT", "Length"]
    outdata.to_csv("SequenceDataDistribution.csv") # save the data distribution of all the sequences

    DataVisualisation("SequenceDataDistribution.csv")
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



def CheckDataDistribution(Input_Folder: str, entry, atmin, atmax, minlength):
    print("----- Processing Data -----")
    ATandLengthDistribution = []

    Total_Accepted_Bases: int = 0
    for Input_File in os.listdir(str(Input_Folder)):
        subprocess.run(["gzip", "-dk", Input_Folder+"/"+Input_File]) # Unzip the raw read files

        with open(Input_Folder+"/"+Input_File[:-3]) as file: # Open the file
            raw_data: list[str] = file.readlines() # Read the file as a list of lines
            location: int = 0 # Keep track of which line in the file I am
            for line in raw_data: # For each line
                if line[0] == "@": # Does it start with an identifier of a new read?
                    sequence: str = raw_data[location+1] # look at the next line, this is the sequence!
                    sequence_length: int = 0 # keep track of total length of sequence
                    at_content: int = 0 # Keep track of the number of Adenines or Thymines encountered
                    for base in sequence:
                        sequence_length += 1
                        if base == "A" or base == "T" or base == "U":
                            at_content += 1
                    at_percentage: float = (at_content/sequence_length) * 100
                    if sequence_length > minlength:
                        ATandLengthDistribution.append((at_percentage,sequence_length))
                        if at_percentage >= atmin and at_percentage < atmax:
                            Total_Accepted_Bases += sequence_length
                location += 1
        os.remove(Input_Folder + "/" + Input_File[:-3]) # Delete the un-compressed file

    if input_settings["genomepath"] != "":
        genomesize: int = 0
        with open(input_settings["genomepath"]) as file:
            for line in file.readlines():
                if line[0] != ">":
                    for character in line:
                        genomesize += 1
        coverage: float = Total_Accepted_Bases / genomesize
        entry.config(state=tk.NORMAL)
        entry.delete('1.0', tk.END)
        entry.insert(tk.END, str(coverage))
        entry.config(state=tk.DISABLED)

    if os.path.isfile("SequenceDataDistribution.csv"):
        os.remove("SequenceDataDistribution.csv")
    outdata = pd.DataFrame(ATandLengthDistribution)
    outdata.columns = ["AT", "Length"]
    outdata.to_csv("SequenceDataDistribution.csv") # save the data distribution of all the sequences
  
def DataVisualisation(input_file): # Plot a nice graph of AT distrubution

    print("----- Creating Visualisation -----")

    data = pd.read_csv(input_file)
    at_percentages: list = data["AT"].to_list()
    sequence_length: list = data["Length"].to_list()

    thresholds = np.arange(0,101,1) # Custom histogram X axis

    ATbins = np.zeros(100)
    Lengthbins = np.zeros(100)

    for index, value in enumerate(at_percentages):
        for t in thresholds: # from 0 to 100 % AT
            if t < 100:
                if (value < t) and (value >= t-1): # if AT% of sequence is less then the upper limit but greater then the lower limit (prior bin)
                    ATbins[t] += 1 # Add one member to this AT% bin
                    Lengthbins[t] += sequence_length[index] # Add up the total number of bases in this bin!
                    break # Get out of loop checking bins, move on to next sequence

    fig, ax = plt.subplots()

    ax.set_xlim([0,100]) # AT% from 0 to 100 %

    ax.axvspan(0, float(ATmin.get()), alpha=0.2, color="red")
    ax.axvspan(float(ATmin.get()), float(ATmax.get()), alpha=0.2, color="green")
    ax.axvline(x=float(ATmin.get()), color="black", ls="--", label="Plastome\nThreshold")

    ax.bar(thresholds[:-1], ATbins*Lengthbins, width=1, color="royalblue")

    ax.set_ylabel("Frequency Density times Length of read", family="serif")
    ax.set_xlabel("AT percentage of read", family="serif")

    plt.legend(prop={"family": "serif"})
    plt.title("AT % distribution of dataset with applied filters", family="serif")
    plt.savefig('DataDistribution.png', bbox_inches='tight')

    img = tk.PhotoImage(file = "DataDistribution.png")
    imageLabel.configure(image = img)
    imageLabel.image = img

    print("----- Finished Data Visualisation -----")

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

def CheckData():
    CheckDataDistribution(input_settings["inputdir"], Coverage, float(ATmin.get()), float(ATmax.get()), float(minReadLength.get()))
    DataVisualisation("SequenceDataDistribution.csv")
    
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


inputframe = tk.Frame(root, borderwidth=1, relief="solid", pady=20, padx=20)
inputframe.grid(column=0, row=1, pady=10, padx=10)
visframe = tk.Frame(root, borderwidth=1, relief="sunken", pady=20, padx=20)
visframe.grid(column=1, row=1, columnspan=2, pady=10, padx=10)
coverageframe = tk.Frame(visframe)
pipelineframe = tk.Frame(root, borderwidth=1, relief="solid", pady=20, padx=20)
pipelineframe.grid(column=0, row=2, columnspan=3, pady=10, padx=10)

TitleLabel = tk.Label(root, text="Plastome Pipeline").grid(row=0, column=1, pady=20)

inputfolderlabel = tk.Label(inputframe, text="Select Input Folder:").grid(row=0, column=0)
SelectInput = tk.Button(inputframe, text="Input Folder", command=OpenInputFolder).grid(row=1, column=0)

minATlabel = tk.Label(inputframe, text="AT % minimum:").grid(row=2, column=0)
MinAT = tk.Entry(inputframe, textvariable=ATmin)
MinAT.insert(0,"58")
MinAT.grid(row=3, column=0)
maxATlabel = tk.Label(inputframe, text="AT % maximum:").grid(row=4, column=0)
MaxAT = tk.Entry(inputframe, textvariable=ATmax)
MaxAT.insert(0,"100")
MaxAT.grid(row=5, column=0)
minreadlengthlabel = tk.Label(inputframe, text="minimum read length (bases):").grid(row=6, column=0)
MinReadLength = tk.Entry(inputframe, textvariable=minReadLength)
MinReadLength.insert(0,"10000")
MinReadLength.grid(row=7, column=0)

referencegenomelabel = tk.Label(inputframe, text="Select reference genome (Align only):").grid(row=8, column=0)
SelectInputGenome = tk.Button(inputframe, text="Reference Genome", command=GenomeLocation).grid(row=9, column=0)

DataVis = tk.Button(visframe, text="Check Data Distribution", command=CheckData).grid(row=0, column=0)
coverageframe.grid(row=1, column=0)
coveragelabel = tk.Label(coverageframe, text="Coverage:").grid(row=0, column=0)
Coverage = tk.Text(coverageframe, height=1, width=10)
Coverage.config(state=tk.DISABLED)
Coverage.grid(row=0, column=1)

startprocesseslabel = tk.Label(pipelineframe, text="Run pipelines:").grid(row=0, column=1,)
LaunchClean = tk.Button(pipelineframe, text="ATCleaner", command=Clean).grid(row=1, column=0, padx=10,pady=5)
CleanCheck = tk.Checkbutton(pipelineframe, text="Has AT cleaner run?", variable=ATcleanerhasrun).grid(row=2, column=0,padx=10)
LaunchAlign = tk.Button(pipelineframe, text="Alignment Pipeline", command=Alignment).grid(row=1, column=1,padx=5,pady=5)
LaunchAssemble = tk.Button(pipelineframe, text="Assembly Pipeline", command=Assembly).grid(row=1, column=2,padx=5,pady=5)

imageLabel = tk.Label(visframe)
imageLabel.grid(row=2, column=0,)

root.mainloop()
