#!/bin/bash

input_folder=$1$"*"
min_length=$2
min_quality=$3
ATmin=$4
ATmax=$5

for file in $input_folder;
do
gzip -dk $file
./data_collection ${file%.gz} $min_length $min_quality > processed_read_data.csv;
rm ${file%.gz}
done

./prepare_histogram processed_read_data.csv > data_to_plot.csv

gnuplot -p -c gnuplot_histogram $ATmin $ATmax data_to_plot.csv

rm processed_read_data.csv
rm data_to_plot.csv
