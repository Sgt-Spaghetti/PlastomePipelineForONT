import sys
import os
import shlex
import shutil
import gzip
import subprocess
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import numpy as np
import statistics
import matplotlib.pyplot as plt

'''
Housekeeping
'''

current_directory: str = os.getcwd()

'''
AT CLEANING ALGORITHM
'''

def ATCLEANER(InputFolder: str, OutputFolder: str, CutoffMin: float, CutoffMax: float, MinRead: float, minquality):

    print("------- STARTING ATCLEANER ---------")

    Input_Folder: str = os.fspath(InputFolder) # Name of the folder holding all the raw reads
    Cutoff_Percentage_From: float  = CutoffMin # Discard entries with an AT percentage below this ammount
    Cutoff_Percentage_To: float  = CutoffMax # Discard entries with an AT percentage above this ammount
    Minimum_Read_Length: float = MinRead # Di truescard entries with a read length below this value
    Minimum_Sequence_Quality: float = minquality
    if OutputFolder == "":
        Output_Folder: str = os.fspath(current_directory + "/ATCleanerOutput") # name of output directory
    else:
        Output_Folder:str = os.fspath(OutputFolder + "/ATCleanerOutput")


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

    Total_Accepted_Bases: int = 0 # how many bases have we got in our data?
    print("------- STARTED CLEANING ---------")

    for Input_File in os.listdir(str(Input_Folder)):
        print("Processing File: " + str(Input_File))
        subprocess.run(["gzip -dk " + (Input_Folder+"/"+Input_File).replace(" ", "\\ ")], shell=True) # Unzip the raw read files
        if UseBinary.get() == True:
            if str(Architecture.get()) == "linux x86_64":
                subprocess.run(["./ATcleaning_linux_x86_64 " + str(Input_Folder+"/"+Input_File[:-3]).replace(" ", "\\ ") +" "+ str(Cutoff_Percentage_From) +" "+ str(Cutoff_Percentage_To) +" "+ str(Minimum_Read_Length) +" "+ str(Minimum_Sequence_Quality) +" "+ str(Output_Folder.replace(" ", "\\ ") + "/" + Input_File.rpartition("/")[2][0:-9] + "_AT-cleaned.fastq")], shell=True)
            elif str(Architecture.get()) == "Mac ARMv8.6A":
                subprocess.run(["./ATcleaning_mac_ARMv8.6A " + str(Input_Folder+"/"+Input_File[:-3]).replace(" ", "\\ ") +" "+ str(Cutoff_Percentage_From) +" "+ str(Cutoff_Percentage_To) +" "+ str(Minimum_Read_Length) +" "+ str(Minimum_Sequence_Quality) +" "+ str(Output_Folder.replace(" ", "\\ ") + "/" + Input_File.rpartition("/")[2][0:-9] + "_AT-cleaned.fastq")], shell=True)
            os.remove(Input_Folder + "/" + Input_File[:-3]) # Delete the un-compressed file
        else:
            with open(Input_Folder +"/"+Input_File[:-3]) as file: # Open the file
                Total_Entries: int = 0 # How many total reads
                Accepted_Entries: int = 0 # How any chloroplast reads have we kept?
                Accepted_Reads: list[str] = []
                current_entry: list[str] = ["","","",""]
                raw_data: list[str] = file.readlines() # Read the file as a list of lines
                location: int = 0 # Keep track of which line in the file I am
                keep_read: bool = True
                lines_per_read: int  = 4

                at_percentage: float = 0.0
                sequence_length: int = 0

                for line in raw_data: # For each line
                    
                    if location < lines_per_read - 1:

                        current_entry[location] = line

                        if location == 1: # location of the sequence
                            sequence: str = line.rstrip()
                            sequence_length = len(sequence)

                            if sequence_length < Minimum_Read_Length:
                                keep_read = False
                            else:
                                keep_read = True
                                at_content: int = 0 # Keep track of the number of Adenines or Thymines encountered
                                for base in sequence:
                                    if base == "A" or base == "T" or base == "U":
                                        at_content += 1
                                at_percentage = (at_content/sequence_length) * 100
                                if at_percentage < Cutoff_Percentage_From or at_percentage >= Cutoff_Percentage_To:
                                    keep_read = False
                                else: keep_read = True
                        location += 1
                    else: # We are at location 4, sequence quality score
                        current_entry[location]: str = line
                        qscore_sequence: str = line.rstrip()
                        if keep_read == True:
                            cumulative_quality_score: int = 0
                            for score in qscore_sequence:
                                cumulative_quality_score += ord(score) - 33 # Phred score is ascii code - 33
                            Qmean: int =  int((cumulative_quality_score / sequence_length)+0.5)
                            if Qmean < Minimum_Sequence_Quality:
                                keep_read = False

                            else:
                                keep_read = True
                                Accepted_Entries += 1 # We have accepted this sequence!
                                for e in current_entry:
                                    Accepted_Reads.append(e)
                                ATandlengthDistribution.append((at_percentage, sequence_length, Qmean, "-")) # Keep track of the length and AT content of all sequences for data visualisation
                                Total_Accepted_Bases += sequence_length
                        
                        Total_Entries += 1 
                        location = 0

                if Accepted_Entries > 0:
                    with open(Output_Folder.replace(" ", "\\ ") + "/" + Input_File.rpartition("/")[2][0:-9] + "_AT-cleaned.fastq", "w") as output_file:
                        for r in Accepted_Reads:
                            output_file.write(r)
                    with open(Output_Folder.replace(" ", "\\ ")+"/ATCLEANER_LOG_"+Input_File.rpartition("/")[2][0:-9] + ".txt", "w") as logfile:
                        logfile.writelines("Total Reads: " + str(Total_Entries) + "\n")
                        logfile.writelines("Accepted Reads: " + str(Accepted_Entries) + "\n")
                        logfile.writelines("Percentage of total kept: " + str((Accepted_Entries/Total_Entries)*100) + "\n")
                        logfile.writelines("Output file name: " + Output_Folder + Input_File[0:-9] + "_AT-cleaned.fastq" + "\n\n")

                    '''
                                with open(Output_Folder+"/ATCLEANER_LOG_"+Input_File.rpartition("/")[2][0:-9] + ".txt", "a") as logfile:
                                    logfile.writelines("entry number " + str(Total_Entries) + " AT percentage: " + str(at_percentage) + " Read Length: " + str(sequence_length) + "\n")
                                    with open(Output_Folder+"/ATCLEANER_LOG_"+Input_File.rpartition("/")[2][0:-9] + ".txt", "a") as logfile: # Create a log file entry for debugging purposes
                                        logfile.writelines("ACCEPTED\n")
                            else: # Otherwise, the constraints were not satisfied
                                with open(Output_Folder+"/ATCLEANER_LOG_"+Input_File.rpartition("/")[2][0:-9] + ".txt", "a") as logfile: # Log the sequence as rejected
                                    logfile.writelines("REJECTED\n")
                    '''

                os.remove(Input_Folder + "/" + Input_File[:-3]) # Delete the un-compressed file

    print("------- FINISHED CLEANING ---------")
    print("------- COMPRESSING OUTPUT---------")
    subprocess.run(["gzip -r " + Output_Folder.replace(" ", "\\ ") + "/*.fastq"], shell=True) # Compress the output files
    ATcleanerhasrun.set(True)

    if UseBinary.get() == False:
        if input_settings["genomepath"] != "":
            genomesize: int = 0
            with open(os.fspath(input_settings["genomepath"])) as file:
                for line in file.readlines():
                    if line[0] != ">":
                        for character in line:
                            genomesize += 1

            coverage: float = Total_Accepted_Bases / genomesize
            Coverage.config(state=tk.NORMAL)
            Coverage.delete('1.0', tk.END)
            Coverage.insert(tk.END, str(coverage))
            Coverage.config(state=tk.DISABLED)
            ATandlengthDistribution.append(("","","",coverage))
                            

        if os.path.isfile(os.path.split(Output_Folder)[0]+"/SequenceDataDistribution.csv"):
            os.remove(os.path.split(Output_Folder)[0]+"/SequenceDataDistribution.csv")
        outdata = pd.DataFrame(ATandlengthDistribution)
        outdata.columns = ["AT", "Length", "Score", "coverage"]
        outdata.to_csv(os.path.split(Output_Folder)[0]+"/SequenceDataDistribution.csv") # save the data distribution of all the sequences

        DataVisualisation(os.path.split(Output_Folder)[0]+"/SequenceDataDistribution.csv", os.path.split(Output_Folder)[0])
    print("------- ATCLEANER FINISHED ---------")


def PlastomeAlignment(ReferenceGenomePath: str, OutputFolder: str):
    ReferenceGenomePath = os.fspath(ReferenceGenomePath)
    if OutputFolder == "":
        Output_Folder: str = os.fspath(current_directory + "/PlastomeAlignmentOutput") # name of output directory
        Input_Folder: str = os.fspath(current_directory + "/ATCleanerOutput") # name of output directory
    else:
        Output_Folder: str = os.fspath(OutputFolder + "/PlastomeAlignmentOutput")
        Input_Folder: str = os.fspath(OutputFolder + "/ATCleanerOutput") # name of output directory

    if os.path.isdir(Output_Folder):
        for f in os.listdir(Output_Folder):
            os.remove(Output_Folder + "/" + f)
        os.removedirs(Output_Folder)

        os.makedirs(Output_Folder)
    else:
        os.makedirs(Output_Folder)

    subprocess.run(["minimap2 -ax map-ont " + ReferenceGenomePath.replace(" ", "\\ ") + " " + Input_Folder.replace(" ", "\\ ") + "/*.fastq.gz" + " > " + Output_Folder.replace(" ", "\\ ") + "/alignment.sam"], shell=True)

    subprocess.run(["samtools", "view", "-b", "-o", Output_Folder+"/alignment.bam", Output_Folder+"/alignment.sam"])

    subprocess.run(["samtools", "sort", "-O", "bam", "-o", Output_Folder+"/sorted_alignment.bam", Output_Folder+"/alignment.bam"])

    subprocess.run(["samtools", "index", Output_Folder+"/sorted_alignment.bam"])

    os.remove(Output_Folder+"/alignment.bam")
    os.remove(Output_Folder+"/alignment.sam")

def PlastomeAssemble(ReferenceGenomePath: str, OutputFolder: str, FlyeParameters: str):
    ReferenceGenomePath = os.fspath(ReferenceGenomePath)
    if OutputFolder == "":
        Output_Folder: str = os.fspath(current_directory + "/PlastomeDeNovoAssemblyOutput")# name of output directory
        Input_Folder: str = os.fspath(current_directory + "/ATCleanerOutput") # name of output directory
    else:
        Output_Folder: str = os.fspath(OutputFolder + "/PlastomeDeNovoAssemblyOutput")
        Input_Folder: str = os.fspath(OutputFolder + "/ATCleanerOutput") # name of output directory

    if os.path.isdir(Output_Folder):
        shutil.rmtree(Output_Folder)
        os.makedirs(Output_Folder)
        os.makedirs(Output_Folder+"/AlignmentCheck")
        os.makedirs(Output_Folder+"/PorechopCleaned")
    else:
        os.makedirs(Output_Folder)
        os.makedirs(Output_Folder+"/AlignmentCheck")
        os.makedirs(Output_Folder+"/PorechopCleaned")

    if UseAlignment.get() == True:
        subprocess.run(["samtools view -bF 0x800" + os.path.split(Output_Folder)[0].replace(" ", "\\ ")+"/PlastomeAlignmentOutput/sorted_alignment.bam > " + Output_Folder.replace(" ", "\\ ") + "/mapped_primaries.bam"], shell=True)
        subprocess.run(["samtools fastq " + Output_Folder.replace(" ", "\\ ") + "/mapped_primaries.bam > " +Output_Folder.replace(" ", "\\ ") +"/mapped_primaries.fastq"], shell=True)
        subprocess.run(["gzip " + Output_Folder.replace(" ", "\\ ") + "/mapped_primaries.fastq"], shell=True)
        if PorechopSkip.get() == False: # If we ARE running porechop, run it on the concatenated aligned reads extracted previously.
            subprocess.run(["porechop -i " + Output_Folder.replace(" ", "\\ ") + "/mapped_primaries.fastq.gz" + " -o " + Output_Folder.replace(" ", "\\ ")+"/concatenated_porechop_cleaned.fastq.gz"], shell=True)
        else: # Else, simply copy the previously extracted aligned reads to a new location
            subprocess.run(["cp "+Output_Folder.replace(" ", "\\ ")+"/mapped_primaries.fastq.gz "+ Output_Folder.replace(" ", "\\ ")+"/concatenated_porechop_cleaned.fastq.gz"], shell=True)
        os.remove(Output_Folder+"/mapped_primaries.bam")
        os.remove(Output_Folder+"/mapped_primaries.fastq.gz")
    else: # If we are not using the prevously aligned reads (default)
        if PorechopSkip.get() == False: # If we are NOT skipping porechop
            for file in os.listdir(Input_Folder):
                if file[-9:] == ".fastq.gz":
                    subprocess.run(["porechop -i " + Input_Folder.replace(" ", "\\ ")+"/"+file + " -o " + Output_Folder.replace(" ", "\\ ")+"/PorechopCleaned/PORECHOP_"+file], shell=True)
                    subprocess.run(["cat " + Output_Folder.replace(" ", "\\ ")+"/PorechopCleaned/*.fastq.gz > " + Output_Folder.replace(" ", "\\ ")+"/concatenated_porechop_cleaned.fastq.gz"], shell=True)
        else:
            subprocess.run(["cat "+Input_Folder.replace(" ", "\\ ")+"/*.fastq.gz > "+ Output_Folder.replace(" ", "\\ ")+"/concatenated_porechop_cleaned.fastq.gz"], shell=True)

    subprocess.run(["flye --nano-hq " + Output_Folder.replace(" ", "\\ ") + "/concatenated_porechop_cleaned.fastq.gz --out-dir " + Output_Folder.replace(" ", "\\ ") + "/FlyeOutput " + FlyeParameters], shell=True)
    shutil.rmtree(Output_Folder+"/PorechopCleaned")
    os.remove(Output_Folder+"/concatenated_porechop_cleaned.fastq.gz")
    
    if ReferenceGenomePath != "":
        subprocess.run(["minimap2 -ax map-ont " + ReferenceGenomePath.replace(" ", "\\ ") + " " + Output_Folder.replace(" ", "\\ ")+ "/FlyeOutput/assembly.fasta" + " > " + Output_Folder.replace(" ", "\\ ") + "/AlignmentCheck/alignment.sam"], shell=True)
        subprocess.run(["samtools", "view", "-b", "-o", Output_Folder+"/AlignmentCheck/alignment.bam", Output_Folder+"/AlignmentCheck/alignment.sam"])

        subprocess.run(["samtools", "sort", "-O", "bam", "-o", Output_Folder+"/AlignmentCheck/sorted_alignment.bam", Output_Folder+"/AlignmentCheck/alignment.bam"])

        subprocess.run(["samtools", "index", Output_Folder+"/AlignmentCheck/sorted_alignment.bam"])

        os.remove(Output_Folder+"/AlignmentCheck/alignment.sam")
        os.remove(Output_Folder+"/AlignmentCheck/alignment.bam")


def CheckDataDistribution(Input_Folder: str, Output_folder:str, entry, atmin, atmax, minlength, minquality):
    Input_Folder = os.fspath(Input_Folder)
    if Output_folder == "":
        OutputFolder:str = current_directory;
    else:
        OutputFolder:str = os.fspath(Output_folder)

    print("----- Processing Data -----")
    ATandLengthDistribution = []

    Total_Accepted_Bases: int = 0
    for Input_File in os.listdir(str(Input_Folder)):
        print("Processing File: " + str(Input_File))
        subprocess.run(["gzip -dk " + Input_Folder.replace(" ", "\\ ") + "/"+Input_File], shell=True) # Unzip the raw read files

        with open(os.fspath(Input_Folder+"/"+Input_File[:-3])) as file: # Open the file
            raw_data: list[str] = file.readlines() # Read the file as a list of lines
            location: int = 0 # Keep track of which line in the file I am
            for line in raw_data: # For each line
                if line[0] == "@" and (raw_data[location+1][0] != "@"): # Does it start with an identifier of a new read?
                    sequence: str = raw_data[location+1] # look at the next line, this is the sequence!
                    sequence_length: int = 0 # keep track of total length of sequence
                    at_content: int = 0 # Keep track of the number of Adenines or Thymines encountered
                    #sequence_score: list[int] = [] # Keep track of the per-base Phred scores.
                    sequence_score: int = 0 # Keep track of the per-base Phred scores.
                    for position, base in enumerate(sequence):
                        sequence_length += 1
                        #sequence_score.append(ord(raw_data[location+3][position]) - 33) # Phred score is ascii code - 33
                        sequence_score += ord(raw_data[location+3][position]) - 33 # Phred score is ascii code - 33
                        if base == "A" or base == "T" or base == "U":
                            at_content += 1
                    at_percentage: float = (at_content/sequence_length) * 100
                    #quality_score: int = statistics.median(sequence_score)
                    quality_score: float = sequence_score/sequence_length
                    if sequence_length > minlength and quality_score >= minquality:
                        ATandLengthDistribution.append((at_percentage,sequence_length,quality_score, "-"))
                        if at_percentage >= atmin and at_percentage < atmax:
                            Total_Accepted_Bases += sequence_length
                location += 1
        os.remove(Input_Folder + "/" + Input_File[:-3]) # Delete the un-compressed file
    if input_settings["genomepath"] != "":
        genomesize: int = 0
        with open((input_settings["genomepath"])) as file:
            for line in file.readlines():
                if line[0] != ">":
                    for character in line:
                        genomesize += 1
        coverage: float = Total_Accepted_Bases / genomesize
        entry.config(state=tk.NORMAL)
        entry.delete('1.0', tk.END)
        entry.insert(tk.END, str(coverage))
        entry.config(state=tk.DISABLED)
        ATandLengthDistribution.append(("", "", "", coverage))

    if os.path.isfile(OutputFolder + "/SequenceDataDistribution.csv"):
        os.remove(OutputFolder + "/SequenceDataDistribution.csv")
    outdata = pd.DataFrame(ATandLengthDistribution)
    outdata.columns = ["AT", "Length", "Score", "coverage"]
    outdata.to_csv(OutputFolder + "/SequenceDataDistribution.csv") # save the data distribution of all the sequences
  
def DataVisualisation(input_file, Output_Folder): # Plot a nice graph of AT distrubution
    if Output_Folder == "":
        OutputFolder:str = current_directory
    else:
        OutputFolder:str = os.fspath(Output_Folder)


    print("----- Creating Visualisation -----")

    data = pd.read_csv(input_file)
    at_percentages: list = data["AT"].to_list()
    sequence_length: list = data["Length"].to_list()

    thresholds = np.arange(0,101,1) # Custom histogram X axis

    #ATbins = np.zeros(100)
    Lengthbins = np.zeros(100)

    for index, value in enumerate(at_percentages):
        for t in thresholds: # from 0 to 100 % AT
            if t < 100:
                if (value < t) and (value >= t-1): # if AT% of sequence is less then the upper limit but greater then the lower limit (prior bin)
                    #ATbins[t] += 1 # Add one member to this AT% bin
                    Lengthbins[t] += sequence_length[index] # Add up the total number of bases in this bin!
                    break # Get out of loop checking bins, move on to next sequence

    fig, ax = plt.subplots()

    ax.set_xlim([0,100]) # AT% from 0 to 100 %

    ax.axvspan(0, float(ATmin.get()), alpha=0.2, color="red")
    ax.axvspan(float(ATmin.get()), float(ATmax.get()), alpha=0.2, color="green")
    ax.axvline(x=float(ATmin.get()), color="black", ls="--", label="Plastome\nThreshold")

    ax.bar(thresholds[:-1], Lengthbins, width=1, color="royalblue")

    ax.set_ylabel("Frequency Density times Length of read", family="serif")
    ax.set_xlabel("AT percentage of read", family="serif")

    plt.legend(prop={"family": "serif"})
    plt.title("AT % distribution of dataset with applied filters", family="serif")
    plt.savefig(OutputFolder + '/DataDistribution.png', bbox_inches='tight')

    img = tk.PhotoImage(file = OutputFolder+"/DataDistribution.png")

    imageLabel.configure(image = img)
    imageLabel.image = img

    print("----- Finished Data Visualisation -----")

'''
GUI STUFF
'''

# Requires chmod u+rx

def Clean():
    ATCLEANER(input_settings["inputdir"], input_settings["outputdir"], float(ATmin.get()), float(ATmax.get()), float(minReadLength.get()), float(minQualityScore.get()))

def Alignment():
    print("------- STARTING ALIGNMENT ---------")
    if ATcleanerhasrun.get() == False:
        Clean()
    PlastomeAlignment(input_settings["genomepath"],input_settings["outputdir"]) 
    print("------- FINISHED ALIGNMENT ---------") 

def Assembly():
    if ATcleanerhasrun.get() == False:
        Clean()
    print("-------- STARTING ASSEMBLY ---------") 
    PlastomeAssemble(input_settings["genomepath"],input_settings["outputdir"], FlyeParams.get())
    print("-------- FINISHED ASSEMBLY ---------") 

def OpenInputFolder():
    InputFolder = filedialog.askdirectory()
    input_settings["inputdir"] = InputFolder

def OpenOutputFolder():
    OutputFolder = filedialog.askdirectory()
    input_settings["outputdir"] = OutputFolder

def GenomeLocation():
    ReferenceGenomeLocation = filedialog.askopenfilename(filetypes=[("FASTA","*.fasta *.fa")])
    input_settings["genomepath"] = ReferenceGenomeLocation

def CheckData():
    if ATcleanerhasrun.get() == False:
        CheckDataDistribution(input_settings["inputdir"], input_settings["outputdir"], Coverage, float(ATmin.get()), float(ATmax.get()), float(minReadLength.get()), float(minQualityScore.get()))
    DataVisualisation(input_settings["outputdir"] + "/SequenceDataDistribution.csv", input_settings["outputdir"])
    
root = tk.Tk()

input_settings: dict = {
    "inputdir" : "",
    "outputdir" : "",
    "genomepath" : "",
}

ATmin = tk.StringVar()
ATmax = tk.StringVar()
minReadLength = tk.StringVar()
minQualityScore = tk.StringVar()
ATcleanerhasrun = tk.BooleanVar()
PorechopSkip = tk.BooleanVar()
UseAlignment = tk.BooleanVar()
UseBinary = tk.BooleanVar()
Architecture = tk.StringVar()
Architecture.set("linux x86_64")
FlyeParams = tk.StringVar()
FlyeParams.set("--genome-size 0.0002g --asm-coverage 50 -i 3")

root.title("Plastome Pipeline")


titleframe = tk.Frame(root, borderwidth=1, relief="flat", pady=20, padx=20)
titleframe.grid(column=0, row=0, pady=10, padx=10, columnspan=3)
inputframe = tk.Frame(root, borderwidth=1, relief="solid", pady=20, padx=20)
inputframe.grid(column=0, row=1, pady=10, padx=10)
inputbuttonframe = tk.Frame(inputframe, borderwidth=1, relief="flat", pady=20, padx=20)
inputbuttonframe.grid(column=0, row=0)
visframe = tk.Frame(root, borderwidth=1, relief="sunken", pady=20, padx=20)
visframe.grid(column=1, row=1, columnspan=2, pady=10, padx=10)
coverageframe = tk.Frame(visframe)
pipelineframe = tk.Frame(root, borderwidth=1, relief="solid", pady=20, padx=20)
pipelineframe.grid(column=0, row=2, columnspan=3, pady=10, padx=10)

TitleLabel = tk.Label(titleframe, text="Plastome Pipeline").grid(row=0, column=0, pady=20)
ArchitectureLabel = tk.Label(titleframe, text="Architecture").grid(row=1, column=0, padx=20)
ArchitectureSelect = tk.OptionMenu(titleframe, Architecture, *["linux x86_64", "Mac ARMv8.6A"]).grid(row=2, column=0)

inputfolderlabel = tk.Label(inputbuttonframe, text="Select Input Folder:").grid(row=0, column=0, padx=10)
SelectInput = tk.Button(inputbuttonframe, text="Input Folder", command=OpenInputFolder).grid(row=1, column=0)
outputfolderlabel = tk.Label(inputbuttonframe, text="Select Output Folder:").grid(row=0, column=1, padx=10)
SelectOutput = tk.Button(inputbuttonframe, text="Output Folder", command=OpenOutputFolder).grid(row=1, column=1)

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
minreadqualityscorelabel = tk.Label(inputframe, text="minimum Q score (mean):").grid(row=8, column=0)
MinReadLength = tk.Entry(inputframe, textvariable=minQualityScore)
MinReadLength.insert(0,"30")
MinReadLength.grid(row=9, column=0)

referencegenomelabel = tk.Label(inputframe, text="Select reference genome (Align only):").grid(row=10, column=0)
SelectInputGenome = tk.Button(inputframe, text="Reference Genome", command=GenomeLocation).grid(row=11, column=0)

DataVis = tk.Button(visframe, text="Check Data Distribution", command=CheckData).grid(row=0, column=0)
coverageframe.grid(row=1, column=0)
coveragelabel = tk.Label(coverageframe, text="Coverage:").grid(row=0, column=0)
Coverage = tk.Text(coverageframe, height=1, width=10)
Coverage.config(state=tk.DISABLED)
Coverage.grid(row=0, column=1)

startprocesseslabel = tk.Label(pipelineframe, text="Run pipelines:").grid(row=0, column=1,)
LaunchCleanBtn = tk.Button(pipelineframe, text="ATCleaner", command=Clean).grid(row=1, column=0, padx=10,pady=5)
CleanCheckBtn = tk.Checkbutton(pipelineframe, text="Has AT cleaner run?", variable=ATcleanerhasrun).grid(row=2, column=0,padx=10)
LaunchAlignBrn = tk.Button(pipelineframe, text="Alignment Pipeline", command=Alignment).grid(row=1, column=1,padx=5,pady=5)
LaunchAssembleBtn = tk.Button(pipelineframe, text="Assembly Pipeline", command=Assembly).grid(row=1, column=2,padx=5,pady=5)
FlyOpts = tk.Entry(pipelineframe, textvariable=FlyeParams).grid(row=2, column=2, padx=5, pady=5)
AlignmentUseCheckBtn = tk.Checkbutton(pipelineframe,text="Use aligned reads", variable=UseAlignment).grid(row=3, column=2, padx=5, pady=5)
SkipPorechopCheckBtn = tk.Checkbutton(pipelineframe,text="Skip Porechop", variable=PorechopSkip).grid(row=4, column=2, padx=5, pady=5)
UseBinariesCheckBtn = tk.Checkbutton(pipelineframe, text="Use C++ binary", variable=UseBinary).grid(row=3, column=0,)

imageLabel = tk.Label(visframe)
imageLabel.grid(row=2, column=0,)

root.mainloop()
