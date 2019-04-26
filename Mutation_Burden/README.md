# Influenza Mutation Bureden Analysis
Steps and scripts used for analyzing mutation rate in scIAV RNA-sequencing. 

The steps are as folows:

0) Quality of alignments and duplication levels were observed throughout procedure with FastQC. 

1) Remove PCR duplicates with BBMap (`dedupe_script.sh`).

2) Align de-duplicated reads with Bowtie2 (`alignment_for_DE.pbs`).

3) Use Samtools mpileup to create quantifications of reads across flu genome (`name_of_file`).

4) Python script to parse mpileup files and generate per base measures (`pileup_parser.py` and `wrapper_for_pileup_parser.sh`)

5) Parse CSVs created in step 4 (`Parse_CSV_LP005.R` and `Parse_CSV_LP009.R`)

6) Visualize mutation rate versus depth of coverage (`chr_distributions.ipynb`)
