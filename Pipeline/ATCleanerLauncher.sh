#!/bin/bash

# Input .fastq files should be placed in a folder named RawData in the same location that this file is in. The ATCleaner_V1.py file should also be in the same location as this file. The fastq files may be decompressed from a .gz form by running "gzip -d *" in the RawData folder.

if [ "$#" -ne 4 ]; then
    echo "Usage: bash ATCleanerLauncher.sh *Input Directory Name* *AT% bottom threshold* *AT% top threshold* *Minimum Read Length (bp)*\n eg. bash ATCleanerLauncher.sh RawData 60 100 10000"

else
if [ -d "ATCleanerOutput" ]; then # If an old ouptut directory exists
	rm -rf "ATCleanerOutput" # delete it
	mkdir "ATCleanerOutput" # start afresh
else
	mkdir "ATCleanerOutput" # else, just create it
fi


files=$(find $PWD/$1 -type f -name "*.fastq") # Read the names of all the files in the RawData folder

for file in $files
do 
	python3 ATCleaner_V2.py $file $2 $3 $4 ATCleanerOutput
done
fi
