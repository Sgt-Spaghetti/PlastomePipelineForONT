#!/bin/bash

if [ -d "DeNovoAssembly" ]; then # If an old ouptut directory exists
	rm -rf "DeNovoAssembly" # delete it
    mkdir -p "DeNovoAssembly/AlignmentCheck" # start afresh
    mkdir -p "DeNovoAssembly/PorechopCleaned"
else
    mkdir -p "DeNovoAssembly/PorechopCleaned" # else, just create it
    mkdir -p "DeNovoAssembly/AlignmentCheck"
fi

files=$(find $PWD/ATCleanerOutput/ -type f -name "*.fastq.gz") # Read the names of all the files in the RawData folder

for file in $files
do 

	porechop -i $file -o DeNovoAssembly/PorechopCleaned/PORECHOP_$(basename "$file")
done

cat DeNovoAssembly/PorechopCleaned/*.fastq.gz > DeNovoAssembly/concatenated_porechop_cleaned.fastq.gz

~/Bioinformatics/Flye/bin/flye --nano-hq DeNovoAssembly/concatenated_porechop_cleaned.fastq.gz --out-dir DeNovoAssembly/FLYOUT --genome-size 0.0002g --no-alt-contigs --scaffold

minimap2 DeNovoAssembly/FLYOUT/assembly.fasta DeNovoAssembly/concatenated_porechop_cleaned.fastq.gz > DeNovoAssembly/approx_mapping_for_polishing.paf

racon DeNovoAssembly/concatenated_porechop_cleaned.fastq.gz DeNovoAssembly/approx_mapping_for_polishing.paf DeNovoAssembly/FLYOUT/assembly.fasta > DeNovoAssembly/consensus.fasta

minimap2 -ax map-ont $PWD/Genomes/csb516_cpv4-simp.fasta $PWD/DeNovoAssembly/consensus.fasta > $PWD/DeNovoAssembly/AlignmentCheck/alignment.sam

samtools view -@ n -Sb -o $PWD/DeNovoAssembly/AlignmentCheck/alignment.bam $PWD/DeNovoAssembly/AlignmentCheck/alignment.sam

samtools sort -O bam -o $PWD/DeNovoAssembly/AlignmentCheck/sorted_alignment.bam $PWD/DeNovoAssembly/AlignmentCheck/alignment.bam

samtools index $PWD/DeNovoAssembly/AlignmentCheck/sorted_alignment.bam


