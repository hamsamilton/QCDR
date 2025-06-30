#!/bin/zsh

bam_dir=$1
ref_bed=$2
out_dir=$3
name=$4
# index_bam_files
touch logfile.txt
for bam in "$bam_dir"/*.bam; do 
    echo "// Indexing file $bam" >> logfile.txt
    samtools index "$bam" &
done

wait
echo "All files have been processed" >> logfile.txt

# run RseQC
for bam in "$bam_dir"/*bam; do
    # make output name
    echo "// Running RseQC on $bam" >> logfile.txt
    python geneBody_coverage.py -r $ref_bed -i $bam -o $bam &
done

# move all text files to a separate folder
mkdir out_dir
mkdir out_dir/coverage_outputs
mv *.txt out_dir/coverage_outputs

# remove extraneous files created by rseqc
rm "$bam_dir"/*.bai
rm "$bam_dir"/*.R
rm "$bam_dir"/*.pdf

# run summary code
python summarizeGBCoutputs.py -txtdir coverage_outputs -out out_dir -name name
