#!/bin/bash

#Input 1 for SE reads, 2 for paired end reads (SUPPORT FOR OPTION 2 NOT YET INCLUDED).
#If using paired end reads, input 2 is the length of the reads
#For example, for a single end read sequence you should call the script with
#sh dedupe_script.sh 1

#Requires bbmap and FLASH (Flash only needed for paired end reads)
export PATH=$PATH:/Users/MacProMatt/Desktop/programs/bbmap
export PATH=$PATH:/Users/MacProMatt/Desktop/programs/FLASH-1.2.11
export PATH=$PATH:/Users/MacProMatt/Desktop/programs/fastx_toolkit

cd /Users/MacProMatt/Desktop/Langlois_Project_002

if [ ! -d "deduped" ]; then
  mkdir deduped
fi

if [[ $1 == 1 ]]; then
  for i in *.fastq; do
    clumpify.sh in=${i} out=deduped/${i}.deduped.fastq dedupe
    rm -rf ${i}
  done
else
  echo "Not a valid input"
fi
