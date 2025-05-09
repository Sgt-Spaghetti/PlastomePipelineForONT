#!/bin/bash

# Input .fastq files should be placed in a folder named RawData in the same location that this file is in. The ATCleaner_V1.py file should also be in the same location as this file. The fastq files may be decompressed from a .gz form by running "gzip -d *" in the RawData folder.

if [ -d "ATCleanerOutput" ]; then # If an old ouptut directory exists
	rm -rf "ATCleanerOutput" # delete it
	mkdir "ATCleanerOutput" # start afresh
else
	mkdir "ATCleanerOutput" # else, just create it
fi


files=$(find $PWD/RawData -type f -name "*.fastq") # Read the names of all the files in the RawData folder

for file in $files
do 
	python3 ATCleaner_V2.py $file 55 ATCleanerOutput
done
