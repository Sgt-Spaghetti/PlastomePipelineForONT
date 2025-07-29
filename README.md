# Plastome Pipeline
Semi automated pipeline for building a De Novo plastome sequence using nanopore sequencing data

Requries:
- Python:
  - TKinter package (for GUI)
  - Numpy
  - Pandas
  - Matplotlib
- Porechop
- Samtools
- Minimap2
- Flye

Unfortunately due to limitations of Flye and Porechop, you **cannot** have spaces in the filepaths leading to the output directory.

## How to start
Simply call the script from the terminal using
```
python3 PlastomePipeline.py
```

## To preview the data distribution of the raw ONT reads with applied filters
1. Select the input directory to where the ONT reads are stored in .fastq.gz format
2. Select the desired output directory for the program
3. Set the minimum and maximium AT content for the reads to keep, as well as the minimum read length and read quality score required for positive selection
4. _OPTIONAL_ Select the reference genome of the species (used to calculate coverage)
5. Press "Check Data Distribution"

A .png image of the resulting graph will also be saved in the designated output directory
   
## To filter raw ONT reads by AT content
1. Select the input directory to where the ONT reads are stored in .fastq.gz format
2. Select the desired output directory for the program
3. Set the minimum and maximium AT content for the reads to keep, as well as the minimum read length and read quality score required for positive selection
4. _OPTIONAL_ Select the reference genome of the species (used to calculate coverage)
5. Press the "ATCleaner" button to process the reads

The resulting filtered reads, along with associated log files are stored in the "ATCleaner" folder within the chosen output directory

## To align filtered ONT reads to a reference genome
1. Select the input directory to where the ONT reads are stored in .fastq.gz format
2. Select the desired output directory for the program
3. Set the minimum and maximium AT content for the reads to keep, as well as the minimum read length and read quality score required for positive selection.
4. **REQUIRED** Select the reference genome of the species
5. **IMPORTANT** if AT Cleaner has run previously as a separate step, and you wish to use the output as is without running it again, make sure the "Has AT Cleaner Run" checkbox is ticked. Otherwise, AT Cleaner will automatically run before alignment
6. Press the "Alignment Pipeline" button

The resulting alignment is stored in the "PlastomeAlignmentOutput" folder within the chosen output directory

## To run a De Novo Assembly
### From scratch / default
1. **IMPORTANT** make sure the "Use aligned reads" and "Skip Porechop" checkboxes are **empty** (unticked)
2. Select the input directory to where the ONT reads are stored in .fastq.gz format
3. Select the desired output directory for the program
4. Set the minimum and maximium AT content for the reads to keep, as well as the minimum read length and read quality score required for positive selection.
5. _OPTIONAL_ Select the reference genome of the species (this can be used to automatically align the de novo assembly with the reference genome after everything has finished)
6. **IMPORTANT** if AT Cleaner has run previously as a separate step, and you wish to use the output as is without running it again, make sure the "Has AT Cleaner Run" checkbox is ticked. Otherwise, AT Cleaner will automatically run before alignment
7. Press the "Assembly Pipeline" button

### Using reads previously aligned to a reference genome with this program
1. **IMPORTANT** make sure the "Use aligned reads" checkbox is **filled** (ticked) and the "Skip Porechop" checkbox is empty (unticked)
2. Select the input directory to where the ONT reads are stored in .fastq.gz format
3. Select the desired output directory for the program
4. Set the minimum and maximium AT content for the reads to keep, as well as the minimum read length and read quality score required for positive selection.
5. _OPTIONAL_ Select the reference genome of the species (this can be used to automatically align the de novo assembly with the reference genome after everything has finished)
6. **IMPORTANT** if AT Cleaner has run previously as a separate step, and you wish to use the output as is without running it again, make sure the "Has AT Cleaner Run" checkbox is ticked. Otherwise, AT Cleaner will automatically run before alignment
7. Press the "Assembly Pipeline" button

### To skip the Porechop adapter trimming step
Make sure the "Skip Porechop" checkbox is **filled** (ticked)

The resulting assembly is stored in the "PlastomeDeNovoAssembly/FLYOUT" folder within the chosen output directory, it is called "assembly.fasta". The optional alignment with a reference genome can be found in the "PlastomeDeNovoAssembly/AlignmentCheck" folder within the chosen output directory
