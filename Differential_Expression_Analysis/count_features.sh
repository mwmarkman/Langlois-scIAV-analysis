#!/bin/bash

export PATH=$PATH:/Users/MacProMatt/Desktop/programs/subread-1.6.3-source/bin

cd /Users/MacProMatt/Desktop/mouse/sorted

if [ ! -d "features_count_all" ]; then
	mkdir features_count_all
fi

featureCounts -a '/Users/MacProMatt/Desktop/mouse/Mus_musculus.GRCm38.95.gtf' -o features_count_all/total_file_seq.count --ignoreDup Sample1-1 Sample2-1 Sample3-1 Sample4-1 Sample4-2 Sample4-3 Sample5-1 Sample5-2 Sample5-3 Sample6-1 Sample6-2 Sample6-3
