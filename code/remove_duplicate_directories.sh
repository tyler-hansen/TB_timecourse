#!/bin/bash

#variables
dups="/mnt/c/Users/tylerhansen/OneDrive - The University of Chicago/5_Postdoc/Computation/projects/TB_timecourse/data/run2_duplicates.txt"
run2_dir="/mnt/c/Users/tylerhansen/OneDrive - The University of Chicago/5_Postdoc/Computation/projects/TB_timecourse/data/kallisto_samples_run2"

# Check if the dups file exists
if [ ! -f "$dups" ]; then
    echo "Error: File not found: $dups"
    exit 1
fi

# Check if the directory exists
if [ ! -d "$run2_dir" ]; then
    echo "Error: Directory not found: $run2_dir"
    exit 1
fi

# Read each line from the file and remove the corresponding directory
while IFS= read -r directory; do
    if [ -d "${run2_dir}/${directory}" ]; then
        echo "Removing directory: $directory"
        rm -r "${run2_dir}/${directory}"
    else
        echo "Directory not found: ${run2_dir}/${directory}"
    fi
done < "$dups"

echo "Removal process completed."