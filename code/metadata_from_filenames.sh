#!/bin/bash

# Specify the directory containing the files
directory="/mnt/c/Users/tylerhansen/OneDrive - The University of Chicago/5_Postdoc/Computation/projects/TB_timecourse/data/kallisto_TB-relevant"

# Specify the output file
tmp="/mnt/c/Users/tylerhansen/OneDrive - The University of Chicago/5_Postdoc/Computation/projects/TB_timecourse/results/metadata_tmp.tsv"
output_file="/mnt/c/Users/tylerhansen/OneDrive - The University of Chicago/5_Postdoc/Computation/projects/TB_timecourse/results/metadata.tsv"

# Navigate to the directory. Exit if directory doesn't exist. 
cd "$directory" || exit

#we want the following structure (kallisto_number  Donor_ID  Timepoint  Infection(TB/NI) MOI Location Run  filename  sample_number), this is different for different sample types:

# run2 MONTREAL NI: (write tsv)
ls -1d *NI_MONTREAL_run2 | awk -F_ -v OFS='\t' '{print $1"_"$2,$3,$4,$5,"MOI_NA",$6,$7, $0}' > "$tmp"

# run2 MONTREAL Mtb: (append tsv)
ls -1d *Mtb*_MONTREAL_run2 | awk -F_ -v OFS='\t' '{print $1"_"$2,$3,$4,$5,"MOI_"$7,$8,$9, $0}' >> "$tmp"

# Non-infected runs 1-3. (no location specified)
ls -1d *NI_run{1..3} | awk -F_ -v OFS='\t' '{print $1"_"$2,$3,$4,$5,"MOI_NA","unspecified",$6, $0}' >> "$tmp"

# Mtb runs 1-3. (no location specified)
ls -1d *5_run{1..3} | awk -F_ -v OFS='\t' '{print $1"_"$2,$3,$4,$5,"MOI_"$7,"unspecified",$8, $0}' >> "$tmp"

# NI run 4. (Chicago)
ls -1d *NI-C*_run4 | awk -F_ '{print $1"-"$2"-"$3"-"$0}' | awk -F- -v OFS='\t' '{print $1"_"$2,$3,$4,$5,"MOI_NA",$6,$7,$8"-"$9"-"$10"-"$11"-"$12}' >> "$tmp"

# Mtb run 4. (Chicago)
ls -1d *5-C*_run4 | awk -F_ '{print $1"-"$2"-"$3"-"$0}' | awk -F- -v OFS='\t' '{print $1"_"$2,$3,$4,$5,"MOI_"$7,$8,$9,$10"-"$11"-"$12"-"$13"-"$14"-"$15"-"$16}' >> "$tmp"

# NI run 4. (no location specified)
ls -1d *NI_run4 | awk -F_ '{print $1"-"$2"-"$3"-"$0}' | awk -F- -v OFS='\t' '{print $1"_"$2,$3,$4,$5,"MOI_NA","unspecified",$6, $7"-"$8"-"$9"-"$10}' >> "$tmp"

# Mtb run 4. (no location specified)
ls -1d *5_run4 | awk -F_ '{print $1"-"$2"-"$3"-"$0}' | awk -F- -v OFS='\t' '{print $1"_"$2,$3,$4,$5,"MOI_"$7,"unspecified",$8, $9"-"$10"-"$11"-"$12"-"$13"-"$14}' >> "$tmp"

#add sample ID
awk -F- -v OFS='\t' '{print $0, "sample_"NR}' "$tmp" > "$output_file"

rm "$tmp"