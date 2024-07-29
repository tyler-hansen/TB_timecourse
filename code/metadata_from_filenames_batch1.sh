#!/bin/bash

# Specify the directory containing the files
directory="/mnt/c/Users/tylerhansen/OneDrive - The University of Chicago/5_Postdoc/Computation/projects/TB_timecourse/data/kallisto_batch1_samples/"

# Specify the output file
output_file="/mnt/c/Users/tylerhansen/OneDrive - The University of Chicago/5_Postdoc/Computation/projects/TB_timecourse/data/metadata_batch1.tsv"

# Navigate to the directory. Exit if directory doesn't exist. 
cd "$directory" || exit

#we want the following structure (Donor_ID  Timepoint  Infection filename), this is different for different sample types:

# Non-infected run 1. (no location specified)
ls -1d *NI_run1 | awk -F_ -v OFS='\t' '{print $3,$4,$5,$0}' > "$output_file"

# Mtb run 1. (no location specified)
ls -1d *5_run1 | awk -F_ -v OFS='\t' '{print $3,$4,$5,$0}' >> "$output_file"