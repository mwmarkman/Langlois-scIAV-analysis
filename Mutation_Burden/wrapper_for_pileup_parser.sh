#!/bin/bash

#Set the directory of the folder in which the pileup files are contained here
#INCLUDE TRAILING SLASH OF FILE PATH
location="/Users/MacProMatt/Desktop/pileup_results_LP005/"

cd ${location}

for i in *pileup; do
    j=${location}${i}
    python3 /Users/MacProMatt/Desktop/scripts/Bioinformatics/pileup_parser.py ${j} ${j}.csv
done
