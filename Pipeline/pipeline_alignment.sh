#!/bin/bash

minimap2 -ax map-ont /home/leo/Documents/Work/Algae2025/Genomes/csb516_cpv4-simp.fasta /home/leo/Documents/Work/Algae2025/ATCleanerOutput/* > /home/leo/Documents/Work/Algae2025/alignment.sam

samtools view -@ n -Sb -o alignment.bam alignment.sam

samtools sort -O bam -o sorted_alignment.bam alignment.bam

samtools index sorted_alignment.bam
