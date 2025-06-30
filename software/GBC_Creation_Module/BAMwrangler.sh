#!/bin/bash

# Create the 'collectedBAMS' directory if it doesn't exist
mkdir -p collectedBAMS

# Find files containing '_chr.bam' in their names and move them to 'collectedBAMS'
find . -type f -name '*.bam' -print0 | while IFS= read -r -d $'\0' file; do
    echo "Moving ${file} to collectedBAMS"
    mv "${file}" collectedBAMS/
done
