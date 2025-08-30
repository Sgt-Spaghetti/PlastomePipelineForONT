#!/bin/bash

# Given two alignments, extract the chromosomal location of all reads that map to one alignment 
# in the second alignment

alignment_to_bait=$1
full_alignment=$2

mapped_reads="$(samtools view -F 0x900 -q 1 $alignment_to_bait | awk '{print $1}')"
samtools view -F 0x900 -q 1 $full_alignment > full_alignment

echo $full_reads

for i in $mapped_reads; do
	 #echo $full_reads | awk -v input=$i '$1 == input && $3 != "*" {print $1,$3,$4}'
	 cat full_alignment | awk -v input=$i '$1 == input && $3 != "*" {print $1,$3,$4}'
done

rm full_alignment
