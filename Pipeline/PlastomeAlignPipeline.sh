#!/bin/bash

if [ -d "AlignementToReferenceGenome" ]; then # If an old ouptut directory exists
	rm -rf "AlignementToReferenceGenome" # delete it
    mkdir  "AlignementToReferenceGenome" # start afresh
else
    mkdir "AlignementToReferenceGenome" #else, just create it

fi

if [ "$#" -ne 4 ]; then
    echo "Usage: bash PlastomePipeline.sh *Input Directory Name* *AT% bottom threshold* *AT% top threshold* *Minimum Read Length (bp)*\n eg. bash PlastomePipeline.sh RawData 60 100 10000"
else
    gzip -d $PWD/$1/*
    bash ATCleanerLauncher.sh $1 $2 $3 $4
    gzip ATCleanerOutput/*.fastq
    
    minimap2 -ax map-ont $PWD/Genomes/csb516_cpv4-simp.fasta $PWD/ATCleanerOutput/* > $PWD/AlignementToReferenceGenome/alignment.sam

    samtools view -@ n -Sb -o AlignementToReferenceGenome/alignment.bam AlignementToReferenceGenome/alignment.sam

    samtools sort -O bam -o AlignementToReferenceGenome/sorted_alignment.bam AlignementToReferenceGenome/alignment.bam

    samtools index AlignementToReferenceGenome/sorted_alignment.bam

fi
