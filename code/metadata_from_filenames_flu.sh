#!/bin/bash

# Specify the directory containing the files
directory="/mnt/c/Users/tylerhansen/OneDrive - The University of Chicago/5_Postdoc/Computation/projects/TB_timecourse/data/kallisto_Flu-only"

# Specify the output file
tmp="/mnt/c/Users/tylerhansen/OneDrive - The University of Chicago/5_Postdoc/Computation/projects/TB_timecourse/data/metadata_flu_tmp.tsv"
output_file="/mnt/c/Users/tylerhansen/OneDrive - The University of Chicago/5_Postdoc/Computation/projects/TB_timecourse/data/metadata_flu.tsv"

# Navigate to the directory. Exit if directory doesn't exist. 
cd "$directory" || exit

#we want the following structure (kallisto_number  Donor_ID  Timepoint  Infection(Flu) MOI Location Run  filename  sample_number), this is different for different sample types:

ls -1d *run{2..3}* | awk -F_ -v OFS='\t' '{print $1"_"$2,$3,$4,$5,"MOI_"$7,"unspecified",$8, $0}' >"$tmp"

#add sample ID
awk -F- -v OFS='\t' '{print $0, "sample_"NR}' "$tmp" > "$output_file"

rm "$tmp"