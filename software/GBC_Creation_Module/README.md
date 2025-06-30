# ezGBC

RseQC provides functions for calculating the gene body coverage of samples. However, this function is extremely computationally intensive. This software is a simple wrapper utility designed for using RseQC for gene body coverage on large groups of .bam files that allows for parallelization on unix compute nodes. It also includes a summary function which normalizes and plots the outputs.

How to use:

This software required Python, R and samtools installed.

This software requires a folder containing .bam files and the location of a .bed file. If you need help collecting your bam files. there is a utility BAMwrangrangler.sh which will find all bam files within your current directory and move them into a new folder.

run main.sh with its 4 main inputs, the directory of BAM files, the bedfile location, where the outputs should be saved, and the project name for the summary plots:

The main outputs within the summary folder are.

a .csv file with all of the gene body coverages for all samples 0-1 normalized.
a folder of .txt files with the .txt outputs for each individual normalization.

The consitutent functions can be run separately if required.


