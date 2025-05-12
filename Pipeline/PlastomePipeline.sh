#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: bash PlastomePipeline.sh *Input Directory Name* *AT% bottom threshold* *AT% top threshold* *Minimum Read Length (bp)*\n eg. bash PlastomePipeline.sh RawData 60 100 10000"
else
    bash ATCleanerLauncher.sh $1 $2 $3 $4
    bash pipeline_alignment.sh
fi
