# Differential Experssion Analysis
The code used here was used to analyze differential experssion in FACS sorted populations of infected lung cells.

The procedure is as follows:

0) Quality of alignments and duplication levels were observed throughout procedure with FastQC. 

1) Remove PCR duplicates with BBMap (`dedupe_script.sh`).

2) Align de-duplicated reads with Bowtie2.

3) Quantify transcripts with SubReads (`count_features.sh`)

4) Perform differential expression analysis and normalization with edgeR (`edgeR_from_FC.R`).

5) Corral DE genes with R (`find_DE.R`).  

6) Gene Ontology performed with DAVID online website. 
