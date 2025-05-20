#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: bash PlastomePipeline.sh *Input Directory Name* *AT% bottom threshold* *AT% top threshold* *Minimum Read Length (bp)* *Reference Genome*\n eg. bash PlastomePipeline.sh RawData 60 100 10000"
else

    gzip -d $1
    bash ATCleanerLauncher.sh $1 $2 $3 $4
    gzip ATCleanerOutput/*.fastq
fi

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

flye --nano-hq DeNovoAssembly/concatenated_porechop_cleaned.fastq.gz --out-dir DeNovoAssembly/FLYOUT --genome-size 0.0002g -i 3 --threads 6

minimap2 -ax map-ont $PWD/Genomes/csb516_cpv4-simp.fasta $PWD/DeNovoAssembly/FLYOUT/assembly.fasta > $PWD/DeNovoAssembly/AlignmentCheck/alignment.sam

samtools view -@ n -Sb -o $PWD/DeNovoAssembly/AlignmentCheck/alignment.bam $PWD/DeNovoAssembly/AlignmentCheck/alignment.sam

samtools sort -O bam -o $PWD/DeNovoAssembly/AlignmentCheck/sorted_alignment.bam $PWD/DeNovoAssembly/AlignmentCheck/alignment.bam

samtools index $PWD/DeNovoAssembly/AlignmentCheck/sorted_alignment.bam

rm $PWD/DeNovoAssembly/AlignmentCheck/alignment.bam
rm $PWD/DeNovoAssembly/AlignmentCheck/alignment.sam


